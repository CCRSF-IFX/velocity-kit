"""Batch correction methods that preserve RNA-velocity layer relationships."""

from __future__ import annotations

import copy
import logging
from pathlib import Path
from typing import Optional, Sequence

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp

from .io import read_velocity_input

logger = logging.getLogger(__name__)


def add_arguments(parser) -> None:
    """Add batch-correction arguments to an argparse parser."""
    parser.add_argument("input_path", help="Input .loom or .h5ad file")
    parser.add_argument("-o", "--output", required=True, help="Corrected .h5ad path")
    parser.add_argument(
        "--method",
        choices=("hansen-combat",),
        default="hansen-combat",
        help="Velocity-aware batch correction method (default: hansen-combat)",
    )
    parser.add_argument(
        "--batch-key",
        required=True,
        help="Observation column defining the batches removed by ComBat",
    )
    parser.add_argument(
        "--preserve-key",
        nargs="+",
        default=None,
        metavar="COLUMN",
        help="Biological covariate column(s) preserved in the ComBat design",
    )
    parser.add_argument(
        "--target-sum",
        type=float,
        default=1e4,
        help="Shared S+U library-size target per cell (default: 10000)",
    )
    parser.add_argument(
        "--min-shared-counts",
        type=int,
        default=20,
        help="Minimum total counts required in both S and U (default: 20)",
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=2000,
        help="Batch-aware highly variable genes retained; use 0 for all (default: 2000)",
    )
    parser.add_argument(
        "--no-preserve-zeros",
        action="store_false",
        dest="preserve_zeros",
        help="Allow ComBat to make originally zero S+U entries nonzero",
    )
    parser.set_defaults(preserve_zeros=True)
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=1,
        help="Increase verbosity (-v for INFO, -vv for DEBUG)",
    )


def _sum_axis(matrix, axis: int) -> np.ndarray:
    return np.asarray(matrix.sum(axis=axis)).ravel()


def _row_scale(matrix, scale: np.ndarray):
    if sp.issparse(matrix):
        return sp.diags(scale).dot(matrix).tocsr()
    return np.asarray(matrix, dtype=np.float64) * scale[:, None]


def _log1p_matrix(matrix):
    if sp.issparse(matrix):
        logged = matrix.copy().astype(np.float64)
        logged.data = np.log1p(logged.data)
        return logged
    return np.log1p(np.asarray(matrix, dtype=np.float64))


def _dense(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        return matrix.toarray().astype(np.float64, copy=False)
    return np.asarray(matrix, dtype=np.float64)


def _validate_design(
    adata,
    batch_key: str,
    preserve_keys: Sequence[str],
) -> np.ndarray:
    if batch_key not in adata.obs:
        raise ValueError(f"Batch key {batch_key!r} was not found in adata.obs")
    if adata.obs[batch_key].isna().any():
        raise ValueError(f"Batch key {batch_key!r} contains missing values")
    if batch_key in preserve_keys:
        raise ValueError("The batch key cannot also be a preserved covariate")
    for key in preserve_keys:
        if key not in adata.obs:
            raise ValueError(f"Preserved covariate {key!r} was not found in adata.obs")
        if adata.obs[key].isna().any():
            raise ValueError(f"Preserved covariate {key!r} contains missing values")

    batches = adata.obs[batch_key].astype(str).to_numpy()
    levels, counts = np.unique(batches, return_counts=True)
    if len(levels) < 2:
        raise ValueError(f"Batch correction requires at least two {batch_key!r} levels")
    small = levels[counts < 2]
    if len(small):
        raise ValueError(
            "Every batch must contain at least two cells; too small: "
            f"{', '.join(small)}"
        )
    return levels


def _input_qc(spliced) -> tuple[np.ndarray, np.ndarray]:
    if sp.issparse(spliced):
        genes = np.asarray((spliced > 0).sum(axis=1)).ravel()
    else:
        genes = np.count_nonzero(np.asarray(spliced) > 0, axis=1)
    counts = _sum_axis(spliced, axis=1)
    return genes.astype(np.int64), counts.astype(np.float64)


def hansen_combat_correct(
    adata,
    batch_key: str,
    preserve_keys: Optional[Sequence[str]] = None,
    target_sum: float = 1e4,
    min_shared_counts: int = 20,
    n_top_genes: Optional[int] = 2000,
    preserve_zeros: bool = True,
):
    """Correct total S+U abundance with ComBat and preserve the spliced ratio."""
    preserve_keys = list(dict.fromkeys(preserve_keys or []))
    batch_levels = _validate_design(adata, batch_key, preserve_keys)
    if target_sum <= 0:
        raise ValueError("target_sum must be positive")
    if min_shared_counts < 0:
        raise ValueError("min_shared_counts cannot be negative")
    if n_top_genes is not None and n_top_genes < 0:
        raise ValueError("n_top_genes cannot be negative")

    spliced = adata.layers["spliced"]
    unspliced = adata.layers["unspliced"]
    input_genes_per_cell, input_counts_per_cell = _input_qc(spliced)
    total = spliced + unspliced
    library_size = _sum_axis(total, axis=1).astype(np.float64)
    if np.any(library_size <= 0):
        raise ValueError("Every cell must have at least one spliced or unspliced count")
    scale = target_sum / library_size

    shared_mask = (_sum_axis(spliced, 0) >= min_shared_counts) & (
        _sum_axis(unspliced, 0) >= min_shared_counts
    )
    if not shared_mask.any():
        raise ValueError("No genes pass the shared spliced/unspliced count filter")

    filtered_genes = np.flatnonzero(shared_mask)
    normalized_total = _row_scale(total[:, filtered_genes], scale)
    log_total = _log1p_matrix(normalized_total)
    filtered = ad.AnnData(
        X=log_total,
        obs=adata.obs.copy(),
        var=adata.var.iloc[filtered_genes].copy(),
    )

    requested_top = None if not n_top_genes else int(n_top_genes)
    if requested_top is not None and requested_top < filtered.n_vars:
        sc.pp.highly_variable_genes(
            filtered,
            n_top_genes=requested_top,
            batch_key=batch_key,
        )
        selected_local = np.flatnonzero(filtered.var["highly_variable"].to_numpy())
    else:
        selected_local = np.arange(filtered.n_vars)
    if not len(selected_local):
        raise ValueError("No genes were selected for batch correction")
    selected_genes = filtered_genes[selected_local]

    spliced_raw = _dense(spliced[:, selected_genes])
    unspliced_raw = _dense(unspliced[:, selected_genes])
    total_raw = spliced_raw + unspliced_raw
    ratio = np.divide(
        spliced_raw,
        total_raw,
        out=np.zeros_like(spliced_raw),
        where=total_raw > 0,
    )
    normalized_selected_total = total_raw * scale[:, None]
    combat_input = ad.AnnData(
        X=np.log1p(normalized_selected_total),
        obs=adata.obs.copy(),
        var=adata.var.iloc[selected_genes].copy(),
    )
    sc.pp.combat(
        combat_input,
        key=batch_key,
        covariates=preserve_keys or None,
    )
    corrected_log_total = np.asarray(combat_input.X, dtype=np.float64)
    if not np.isfinite(corrected_log_total).all():
        raise ValueError("ComBat produced non-finite corrected values")

    corrected_total = np.expm1(corrected_log_total)
    zero_mask = total_raw == 0
    negative_mask = corrected_total < 0
    negative_values = int(np.count_nonzero(negative_mask))
    negative_observed = int(np.count_nonzero(negative_mask & ~zero_mask))
    negative_structural_zeros = negative_values - negative_observed
    np.maximum(corrected_total, 0, out=corrected_total)
    if preserve_zeros:
        corrected_total[zero_mask] = 0

    corrected_spliced = corrected_total * ratio
    corrected_unspliced = corrected_total - corrected_spliced
    positive = corrected_total > 0
    ratio_error = 0.0
    if positive.any():
        output_ratio = np.divide(
            corrected_spliced,
            corrected_total,
            out=np.zeros_like(corrected_spliced),
            where=positive,
        )
        ratio_error = float(np.max(np.abs(output_ratio[positive] - ratio[positive])))

    output = ad.AnnData(
        X=np.log1p(corrected_spliced).astype(np.float32),
        obs=adata.obs.copy(),
        var=adata.var.iloc[selected_genes].copy(),
    )
    output.layers["spliced"] = corrected_spliced.astype(np.float32)
    output.layers["unspliced"] = corrected_unspliced.astype(np.float32)
    output.var["highly_variable"] = True
    output.obs["n_genes_by_counts"] = input_genes_per_cell
    output.obs["total_counts"] = input_counts_per_cell
    output.uns["log1p"] = {"base": None}

    for key, value in adata.uns.items():
        if (
            str(key).startswith("velocitykit_")
            and key != "velocitykit_batch_correction"
        ):
            output.uns[key] = copy.deepcopy(value)
    output.uns["velocitykit_batch_correction"] = {
        "method": "hansen-combat",
        "batch_key": batch_key,
        "batch_levels": np.asarray(batch_levels, dtype=str),
        "preserve_keys": np.asarray(preserve_keys, dtype=str),
        "target_sum": float(target_sum),
        "min_shared_counts": int(min_shared_counts),
        "n_top_genes": int(len(selected_genes)),
        "n_input_genes": int(adata.n_vars),
        "n_genes_after_shared_filter": int(len(filtered_genes)),
        "preserve_zeros": bool(preserve_zeros),
        "negative_totals_clipped": negative_values,
        "negative_observed_totals_clipped": negative_observed,
        "negative_structural_zero_totals_clipped": negative_structural_zeros,
        "max_spliced_ratio_error": ratio_error,
    }
    return output


def run(args) -> None:
    """Run Hansen-ComBat correction and write an H5AD file."""
    output_path = Path(args.output)
    if output_path.suffix.lower() != ".h5ad":
        raise ValueError("Batch-corrected output must use the .h5ad extension")

    logger.info("Reading velocity input: %s", args.input_path)
    adata = read_velocity_input(args.input_path)
    corrected = hansen_combat_correct(
        adata,
        batch_key=args.batch_key,
        preserve_keys=args.preserve_key,
        target_sum=args.target_sum,
        min_shared_counts=args.min_shared_counts,
        n_top_genes=args.n_top_genes,
        preserve_zeros=args.preserve_zeros,
    )
    corrected.uns["velocitykit_batch_correction"]["input_path"] = str(
        Path(args.input_path).resolve()
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    corrected.write_h5ad(output_path)
    logger.info(
        "Wrote corrected H5AD with %s cells and %s genes: %s",
        f"{corrected.n_obs:,}",
        f"{corrected.n_vars:,}",
        output_path,
    )
