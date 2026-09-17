"""Assemble multiple velocity inputs into one canonical H5AD file."""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import List, Tuple

import anndata as ad
import numpy as np
import pandas as pd

from .io import read_velocity_input


logger = logging.getLogger(__name__)


def add_arguments(parser) -> None:
    """Add manifest-driven assembly arguments to an argparse parser."""
    parser.add_argument(
        "--manifest",
        required=True,
        help=(
            "CSV or TSV with a required 'path' column and optional per-input "
            "annotations such as sample, batch, and source."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output .h5ad path",
    )
    parser.add_argument(
        "--gene-join",
        choices=("exact", "intersection"),
        default="exact",
        help=(
            "Gene alignment policy: require identical gene sets (exact, default) "
            "or retain their intersection."
        ),
    )
    parser.add_argument(
        "--source-key",
        default="source_file",
        help="Observation column identifying each input (default: source_file)",
    )
    parser.add_argument(
        "--index-separator",
        default=":",
        help="Separator used to append source labels to cell IDs (default: ':')",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=1,
        help="Increase verbosity (-v for INFO, -vv for DEBUG)",
    )


def _read_manifest(manifest_path: str) -> Tuple[pd.DataFrame, Path]:
    path = Path(manifest_path)
    if not path.is_file():
        raise FileNotFoundError(f"Assembly manifest not found: {manifest_path}")
    separator = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    manifest = pd.read_csv(path, sep=separator)
    if manifest.empty:
        raise ValueError("Assembly manifest contains no inputs")
    if "path" not in manifest:
        raise ValueError("Assembly manifest must contain a 'path' column")
    if manifest["path"].isna().any():
        raise ValueError("Assembly manifest 'path' column contains missing values")
    return manifest, path.resolve().parent


def _safe_source_label(value: str) -> str:
    label = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    if not label:
        raise ValueError(f"Invalid empty source label derived from {value!r}")
    return label


def _source_labels(manifest: pd.DataFrame, paths: List[Path]) -> List[str]:
    if "source" in manifest:
        if manifest["source"].isna().any():
            raise ValueError("Assembly manifest 'source' column contains missing values")
        labels = [_safe_source_label(value) for value in manifest["source"]]
    else:
        stems = [_safe_source_label(path.stem) for path in paths]
        counts = pd.Series(stems).value_counts()
        labels = [
            stem if counts[stem] == 1 else f"{stem}_{index + 1}"
            for index, stem in enumerate(stems)
        ]
    if len(set(labels)) != len(labels):
        raise ValueError(
            "Assembly source labels must be unique; provide a unique 'source' column"
        )
    return labels


def _annotate_input(adata, row: pd.Series, source_key: str) -> None:
    for column, value in row.items():
        if column in {"path", "source", source_key}:
            continue
        if pd.isna(value):
            raise ValueError(
                f"Manifest annotation {column!r} contains a missing value"
            )
        if column in adata.obs:
            existing = adata.obs[column].dropna().astype(str).unique()
            if len(existing) and set(existing) != {str(value)}:
                raise ValueError(
                    f"Manifest value {value!r} for {column!r} conflicts with "
                    f"existing input annotations: {', '.join(existing[:5])}"
                )
        adata.obs[column] = value


def assemble_from_manifest(
    manifest_path: str,
    gene_join: str = "exact",
    source_key: str = "source_file",
    index_separator: str = ":",
):
    """Assemble manifest inputs into one canonical in-memory AnnData object."""
    if gene_join not in {"exact", "intersection"}:
        raise ValueError("gene_join must be 'exact' or 'intersection'")
    if not source_key:
        raise ValueError("source_key cannot be empty")
    if not index_separator:
        raise ValueError("index_separator cannot be empty")

    manifest, manifest_dir = _read_manifest(manifest_path)
    paths = []
    for value in manifest["path"].astype(str):
        path = Path(value)
        if not path.is_absolute():
            path = manifest_dir / path
        paths.append(path.resolve())
    labels = _source_labels(manifest, paths)

    inputs = []
    input_n_obs = []
    input_n_vars = []
    for (_, row), path, label in zip(manifest.iterrows(), paths, labels):
        logger.info("Reading assembly input %s: %s", label, path)
        adata = read_velocity_input(str(path))
        if not adata.obs_names.is_unique:
            raise ValueError(f"Input {path} contains duplicate cell identifiers")
        if not adata.var_names.is_unique:
            raise ValueError(f"Input {path} contains duplicate gene identifiers")
        _annotate_input(adata, row, source_key)
        adata.X = adata.layers["spliced"].copy()
        inputs.append(adata)
        input_n_obs.append(int(adata.n_obs))
        input_n_vars.append(int(adata.n_vars))

    reference_genes = list(inputs[0].var_names.astype(str))
    if gene_join == "exact":
        reference_set = set(reference_genes)
        for path, adata in zip(paths[1:], inputs[1:]):
            genes = set(adata.var_names.astype(str))
            if genes != reference_set:
                missing = len(reference_set - genes)
                extra = len(genes - reference_set)
                raise ValueError(
                    f"Gene set for {path} differs from the first input "
                    f"({missing} missing, {extra} extra); use --gene-join intersection"
                )
        genes_to_keep = reference_genes
    else:
        common = set(reference_genes)
        for adata in inputs[1:]:
            common.intersection_update(adata.var_names.astype(str))
        genes_to_keep = [gene for gene in reference_genes if gene in common]
        if not genes_to_keep:
            raise ValueError("Assembly inputs have no genes in common")

    aligned = [adata[:, genes_to_keep].copy() for adata in inputs]
    combined = ad.concat(
        aligned,
        axis=0,
        join="inner",
        merge="same",
        label=source_key,
        keys=labels,
        index_unique=index_separator,
    )
    combined.X = combined.layers["spliced"].copy()
    combined.uns["velocitykit_assembly"] = {
        "manifest": str(Path(manifest_path).resolve()),
        "gene_join": gene_join,
        "source_key": source_key,
        "index_separator": index_separator,
        "input_paths": np.asarray([str(path) for path in paths], dtype=str),
        "source_labels": np.asarray(labels, dtype=str),
        "input_n_obs": np.asarray(input_n_obs, dtype=np.int64),
        "input_n_vars": np.asarray(input_n_vars, dtype=np.int64),
    }
    return combined


def run(args) -> None:
    """Run manifest-driven assembly and write the canonical H5AD output."""
    output = Path(args.output)
    if output.suffix.lower() != ".h5ad":
        raise ValueError("Assembly output must use the .h5ad extension")
    combined = assemble_from_manifest(
        manifest_path=args.manifest,
        gene_join=args.gene_join,
        source_key=args.source_key,
        index_separator=args.index_separator,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(output)
    logger.info(
        "Wrote assembled H5AD with %s cells and %s genes: %s",
        f"{combined.n_obs:,}",
        f"{combined.n_vars:,}",
        output,
    )
