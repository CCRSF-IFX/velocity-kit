"""Generate scVelo analysis reports from loom or H5AD files."""

from __future__ import annotations

import os
import re
import warnings
from datetime import datetime
from html import escape
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt

from .io import read_velocity_input

# Suppress known deprecation warnings from dependencies
warnings.filterwarnings('ignore', category=UserWarning, module='louvain')
warnings.filterwarnings('ignore', message='pkg_resources is deprecated')


def _read_cell_metadata(metadata_file: str) -> pd.DataFrame:
    """Read comma- or tab-delimited cell metadata."""
    path = Path(metadata_file)
    if not path.is_file():
        raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

    separator = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    metadata = pd.read_csv(path, sep=separator)
    if metadata.empty:
        raise ValueError(f"Metadata file contains no rows: {metadata_file}")
    return metadata


def _values_match(left: pd.Series, right: pd.Series) -> bool:
    """Return whether two aligned metadata columns agree on non-missing values."""
    compared = left.notna() & right.notna()
    if not compared.any():
        return True

    left_values = left[compared]
    right_values = right[compared]
    if pd.api.types.is_numeric_dtype(left_values) and pd.api.types.is_numeric_dtype(
        right_values
    ):
        return bool(
            np.allclose(
                left_values.to_numpy(dtype=float),
                right_values.to_numpy(dtype=float),
                equal_nan=True,
            )
        )
    return bool(
        np.array_equal(
            left_values.astype(str).to_numpy(), right_values.astype(str).to_numpy()
        )
    )


def _adata_key_values(adata, key: str) -> pd.Series:
    """Get a cell identifier from ``adata.obs`` or its index."""
    if key == "_index":
        return pd.Series(adata.obs_names.astype(str), index=adata.obs_names)
    if key not in adata.obs:
        raise ValueError(
            f"AnnData cell key {key!r} was not found in adata.obs. "
            "Use '_index' to join against adata.obs_names."
        )
    values = adata.obs[key]
    if values.isna().any():
        raise ValueError(f"AnnData cell key {key!r} contains missing values")
    return pd.Series(values.astype(str).to_numpy(), index=adata.obs_names)


def _resolve_metadata_keys(
    adata,
    metadata: pd.DataFrame,
    metadata_key: Optional[str],
    adata_key: Optional[str],
) -> Tuple[str, str]:
    """Resolve a unique metadata-to-AnnData cell identifier mapping."""
    if adata_key is not None and metadata_key is None:
        if adata_key == "_index":
            raise ValueError("--metadata-key is required with --adata-key _index")
        metadata_key = adata_key

    if metadata_key is not None:
        if metadata_key not in metadata:
            raise ValueError(
                f"Metadata key {metadata_key!r} was not found. Available columns: "
                f"{', '.join(map(str, metadata.columns))}"
            )
        if adata_key is None:
            if metadata_key in adata.obs:
                adata_key = metadata_key
            else:
                metadata_values = set(metadata[metadata_key].dropna().astype(str))
                if set(adata.obs_names.astype(str)).issubset(metadata_values):
                    adata_key = "_index"
                else:
                    raise ValueError(
                        f"Could not match metadata key {metadata_key!r} to adata.obs. "
                        "Provide --adata-key explicitly."
                    )
        return metadata_key, adata_key

    candidates = []
    for column in metadata.columns:
        if column not in adata.obs:
            continue
        metadata_values = metadata[column]
        adata_values = adata.obs[column]
        if (
            metadata_values.notna().all()
            and adata_values.notna().all()
            and metadata_values.astype(str).is_unique
            and adata_values.astype(str).is_unique
            and set(adata_values.astype(str)).issubset(
                set(metadata_values.astype(str))
            )
        ):
            candidates.append((str(column), str(column)))

    for column in metadata.columns:
        metadata_values = metadata[column]
        if (
            metadata_values.notna().all()
            and metadata_values.astype(str).is_unique
            and set(adata.obs_names.astype(str)).issubset(
                set(metadata_values.astype(str))
            )
        ):
            candidates.append((str(column), "_index"))

    candidates = list(dict.fromkeys(candidates))
    if len(candidates) != 1:
        candidate_text = ", ".join(f"{m}:{a}" for m, a in candidates) or "none"
        raise ValueError(
            "Could not infer one unique metadata join. Provide --metadata-key and "
            f"--adata-key explicitly. Candidate mappings: {candidate_text}"
        )
    return candidates[0]


def attach_cell_metadata(
    adata,
    metadata_file: str,
    metadata_key: Optional[str] = None,
    adata_key: Optional[str] = None,
) -> Dict[str, Any]:
    """Attach external cell metadata after a validated one-to-one join."""
    metadata = _read_cell_metadata(metadata_file)
    metadata_key, adata_key = _resolve_metadata_keys(
        adata, metadata, metadata_key, adata_key
    )

    metadata_ids = metadata[metadata_key]
    if metadata_ids.isna().any():
        raise ValueError(f"Metadata key {metadata_key!r} contains missing values")
    metadata_ids = metadata_ids.astype(str)
    if not metadata_ids.is_unique:
        examples = metadata_ids[metadata_ids.duplicated(keep=False)].unique()[:5]
        raise ValueError(
            f"Metadata key {metadata_key!r} contains duplicate cell identifiers: "
            f"{', '.join(examples)}"
        )

    adata_ids = _adata_key_values(adata, adata_key)
    if adata_ids.isna().any() or not adata_ids.is_unique:
        raise ValueError(f"AnnData cell key {adata_key!r} must be complete and unique")

    indexed_metadata = metadata.copy()
    indexed_metadata.index = metadata_ids
    missing = adata_ids[~adata_ids.isin(indexed_metadata.index)]
    if len(missing):
        raise ValueError(
            f"Metadata is missing {len(missing):,} of {adata.n_obs:,} cells; "
            f"examples: {', '.join(missing.iloc[:5])}"
        )

    aligned = indexed_metadata.loc[adata_ids.to_numpy()].copy()
    aligned.index = adata.obs_names
    for column in metadata.columns:
        if column == metadata_key:
            continue
        incoming = aligned[column]
        if column in adata.obs and not _values_match(adata.obs[column], incoming):
            raise ValueError(
                f"External metadata column {column!r} conflicts with existing "
                "adata.obs values after joining cells"
            )
        if pd.api.types.is_object_dtype(incoming) or isinstance(
            incoming.dtype, pd.StringDtype
        ):
            adata.obs[column] = pd.Categorical(
                incoming, categories=pd.unique(incoming.dropna()), ordered=False
            )
        else:
            adata.obs[column] = incoming.to_numpy()

    details = {
        "file": str(metadata_file),
        "metadata_key": metadata_key,
        "adata_key": adata_key,
        "matched_cells": int(adata.n_obs),
        "extra_metadata_rows": int(len(metadata) - adata.n_obs),
        "columns": [str(column) for column in metadata.columns if column != metadata_key],
    }
    adata.uns["velocitykit_external_metadata"] = details
    return details


def _plot_filename(column: str) -> str:
    """Create a safe, stable UMAP filename from an annotation name."""
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", column).strip("._") or "annotation"
    return f"umap_{slug}.png"


def subset_cells(
    adata,
    subset_by: str,
    subset_values: Sequence[str],
) -> Tuple[Any, Dict[str, Any]]:
    """Return cells matching selected observation values and selection details."""
    if subset_by not in adata.obs:
        raise ValueError(
            f"Subset column {subset_by!r} was not found after metadata loading. "
            f"Available columns: {', '.join(map(str, adata.obs.columns))}"
        )

    requested = list(dict.fromkeys(str(value) for value in subset_values))
    if not requested:
        raise ValueError("At least one --subset-values entry is required")

    observed = adata.obs[subset_by]
    available = set(observed.dropna().astype(str))
    missing = [value for value in requested if value not in available]
    if missing:
        available_text = ", ".join(sorted(available))
        raise ValueError(
            f"Requested values were not found in {subset_by!r}: {', '.join(missing)}. "
            f"Available values: {available_text}"
        )

    mask = observed.astype(str).isin(requested).to_numpy()
    original_cells = int(adata.n_obs)
    selected = adata[mask].copy()
    details = {
        "column": str(subset_by),
        "values": requested,
        "original_cells": original_cells,
        "retained_cells": int(selected.n_obs),
    }
    selected.uns["velocitykit_cell_subset"] = details
    return selected, details


def calculate_input_qc_metrics(adata) -> None:
    """Calculate per-cell QC from the input matrix before gene filtering."""
    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
    adata.uns["velocitykit_input_qc"] = {
        "matrix": "X",
        "n_genes_before_filtering": int(adata.n_vars),
    }


def write_analyzed_adata(adata, output_path: str) -> str:
    """Write a fully analyzed AnnData object to an H5AD path."""
    path = Path(output_path)
    if path.suffix.lower() != ".h5ad":
        raise ValueError("Analyzed AnnData output must use the .h5ad extension")
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
    return str(path)


def run_scvelo_and_generate_report(
    input_path: str,
    output_dir: str,
    sample_name: Optional[str] = None,
    metadata_file: Optional[str] = None,
    metadata_key: Optional[str] = None,
    adata_key: Optional[str] = None,
    color_by: Optional[Sequence[str]] = None,
    subset_by: Optional[str] = None,
    subset_values: Optional[Sequence[str]] = None,
    save_anndata: Optional[str] = None,
) -> str:
    """
    Run a standard scVelo pipeline on a loom or H5AD file and generate an HTML report
    with QC and analysis plots.

    Parameters
    ----------
    input_path : str
        Path to an input ``.loom`` or ``.h5ad`` file.
    output_dir : str
        Directory where plots and HTML report will be written.
    sample_name : str, optional
        Name used in plot titles and report filename. If None, derived from input_path.
    metadata_file : str, optional
        CSV or TSV containing cell-level annotations to attach before analysis.
    metadata_key : str, optional
        Unique cell identifier column in the external metadata.
    adata_key : str, optional
        Matching identifier in ``adata.obs``. Use ``_index`` for ``obs_names``.
    color_by : sequence of str, optional
        Annotation columns for additional UMAP plots.
    subset_by : str, optional
        Observation or metadata column used to select cells before preprocessing.
    subset_values : sequence of str, optional
        Values retained from ``subset_by``.
    save_anndata : str, optional
        Output path for the fully analyzed AnnData object.

    Returns
    -------
    report_path : str
        Path to the generated HTML report.
    """

    # -------------------------------------------------------------------------
    # Setup
    # -------------------------------------------------------------------------
    os.makedirs(output_dir, exist_ok=True)
    if sample_name is None:
        sample_name = os.path.splitext(os.path.basename(input_path))[0]

    # make scvelo write figures into output_dir
    scv.settings.figdir = output_dir
    scv.settings.set_figure_params(frameon=False, dpi=120)

    # -------------------------------------------------------------------------
    # Load data
    # -------------------------------------------------------------------------
    print(f"[{sample_name}] Reading velocity input: {input_path}")
    adata = read_velocity_input(input_path)

    color_by = list(dict.fromkeys(color_by or []))
    metadata_details = None
    if metadata_file is not None:
        print(f"[{sample_name}] Attaching cell metadata: {metadata_file}")
        metadata_details = attach_cell_metadata(
            adata,
            metadata_file=metadata_file,
            metadata_key=metadata_key,
            adata_key=adata_key,
        )
        print(
            f"[{sample_name}] Matched {metadata_details['matched_cells']} cells "
            f"using {metadata_details['metadata_key']} -> "
            f"{metadata_details['adata_key']}"
        )

    subset_details = None
    if subset_by is not None:
        adata, subset_details = subset_cells(
            adata,
            subset_by=subset_by,
            subset_values=subset_values or [],
        )
        print(
            f"[{sample_name}] Retained {subset_details['retained_cells']} of "
            f"{subset_details['original_cells']} cells where "
            f"{subset_details['column']} is one of "
            f"{', '.join(subset_details['values'])}"
        )

    missing_colors = [column for column in color_by if column not in adata.obs]
    if missing_colors:
        raise ValueError(
            "Requested --color-by columns were not found after metadata loading: "
            f"{', '.join(missing_colors)}. Available columns: "
            f"{', '.join(map(str, adata.obs.columns))}"
        )
    
    # Fix categorical columns that may cause issues with newer pandas
    # Convert any categorical columns to regular strings to avoid
    # "property 'categories' of 'Categorical' object has no setter" errors
    # for col in adata.obs.columns:
    #     if hasattr(adata.obs[col], 'cat'):
    #         adata.obs[col] = adata.obs[col].astype(str)
    # for col in adata.var.columns:
    #     if hasattr(adata.var[col], 'cat'):
    #         adata.var[col] = adata.var[col].astype(str)

    print(f"[{sample_name}] Data loaded: {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"[{sample_name}] Computing input QC metrics before gene filtering")
    calculate_input_qc_metrics(adata)

    # -------------------------------------------------------------------------
    # Adaptive parameters based on cell count
    # -------------------------------------------------------------------------
    n_cells = adata.n_obs
    
    # Set n_neighbors adaptively (never more than n_cells - 1)
    # Standard: 30, but scale down for small datasets
    if n_cells < 50:
        n_neighbors = max(5, min(15, n_cells - 1))
        print(f"[{sample_name}] Small dataset detected ({n_cells} cells), using n_neighbors={n_neighbors}")
    else:
        n_neighbors = 30
    
    # Set n_pcs adaptively (never more than min(n_cells, n_genes) - 1)
    n_pcs = min(30, n_cells - 1, adata.n_vars - 1)
    if n_pcs < 30:
        print(f"[{sample_name}] Using n_pcs={n_pcs} (limited by dataset size)")

    # -------------------------------------------------------------------------
    # scVelo preprocessing - following scvelo_example.py
    # -------------------------------------------------------------------------
    print(f"[{sample_name}] Running preprocessing and filtering")
    
    # Filter and normalize (combined step from scVelo)
    # scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    # If your data are raw counts in adata.X, this is typical:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)  # <-- replaces scVelo's deprecated log1p

    # Optional but recommended for velocity workflows:
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var["highly_variable"]].copy()

    # PCA + neighbors must be computed explicitly now (scVelo >= 0.4)
    sc.pp.scale(adata, max_value=10)        # optional; many people do it
    # Compute PCA and neighbors (use adaptive parameters)
    sc.tl.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # Compute moments for velocity estimation
    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)

    # -------------------------------------------------------------------------
    # Velocity computation - following scvelo_example.py
    # -------------------------------------------------------------------------
    print(f"[{sample_name}] Computing RNA velocity")
    
    # Run velocity estimation
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)

    # -------------------------------------------------------------------------
    # Embedding (UMAP) - following scvelo_example.py
    # -------------------------------------------------------------------------
    print(f"[{sample_name}] Computing embeddings")
    
    sc.tl.umap(adata)

    # -------------------------------------------------------------------------
    # Optional: Velocity metrics for QC
    # -------------------------------------------------------------------------
    print(f"[{sample_name}] Computing velocity QC metrics")

    # Velocity confidence (optional but useful for QC)
    scv.tl.velocity_confidence(adata)

    # -------------------------------------------------------------------------
    # Plot helpers
    # -------------------------------------------------------------------------
    def save_current_fig(filename: str):
        """Save current matplotlib figure to output_dir/filename and close."""
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, bbox_inches="tight")
        plt.close()
        return filename  # return relative name for HTML

    generated_plots = []

    # -------------------------------------------------------------------------
    # QC plots
    # -------------------------------------------------------------------------
    print(f"[{sample_name}] Generating QC plots")

    # 1. Total counts per cell
    sc.pl.violin(
        adata,
        ["total_counts"],
        jitter=0.4,
        multi_panel=False,
        show=False
    )
    generated_plots.append(
        ("qc_total_counts.png", "QC: Total counts per cell")
    )
    save_current_fig("qc_total_counts.png")

    # 2. Number of genes per cell in the input matrix, before HVG filtering.
    sc.pl.violin(
        adata,
        ["n_genes_by_counts"],
        jitter=0.4,
        multi_panel=False,
        show=False
    )
    generated_plots.append(
        (
            "qc_n_genes_by_counts.png",
            "QC: Number of genes per cell before gene filtering",
        )
    )
    save_current_fig("qc_n_genes_by_counts.png")

    # 3. Spliced/unspliced proportions (if present)
    if {"spliced", "unspliced"}.issubset(adata.layers.keys()):
        scv.pl.proportions(adata, show=False)
        generated_plots.append(
            ("qc_spliced_unspliced_proportions.png",
             "QC: Spliced vs unspliced proportions")
        )
        save_current_fig("qc_spliced_unspliced_proportions.png")

    # -------------------------------------------------------------------------
    # Embedding + velocity plots
    # -------------------------------------------------------------------------
    print(f"[{sample_name}] Generating velocity and embedding plots")

    # 4. UMAP colored by total_counts
    scv.pl.scatter(
        adata,
        basis="umap",
        color="total_counts",
        show=False
    )
    generated_plots.append(
        ("umap_total_counts.png", "UMAP colored by total counts")
    )
    save_current_fig("umap_total_counts.png")

    # 5. UMAP colored by velocity_length
    if "velocity_length" in adata.obs:
        scv.pl.scatter(
            adata,
            basis="umap",
            color="velocity_length",
            show=False
        )
        generated_plots.append(
            ("umap_velocity_length.png", "UMAP colored by velocity length")
        )
        save_current_fig("umap_velocity_length.png")

    # 6. Velocity embedding (arrows on UMAP)
    scv.pl.velocity_embedding(
        adata,
        basis="umap",
        arrow_size=2,
        dpi=120,
        show=False
    )
    generated_plots.append(
        ("velocity_embedding_umap.png",
         "Velocity embedding (arrows) on UMAP")
    )
    save_current_fig("velocity_embedding_umap.png")

    # 7. Velocity stream plot on UMAP
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        dpi=120,
        show=False
    )
    generated_plots.append(
        ("velocity_stream_umap.png",
         "Velocity stream plot on UMAP")
    )
    save_current_fig("velocity_stream_umap.png")

    # Additional UMAPs colored by user-selected cell metadata.
    for column in color_by:
        print(f"[{sample_name}] Plotting UMAP colored by {column}")
        sc.pl.umap(
            adata,
            color=column,
            title=f"UMAP colored by {column}",
            show=False,
        )
        filename = _plot_filename(column)
        generated_plots.append((filename, f"UMAP colored by {column}"))
        save_current_fig(filename)

    # -------------------------------------------------------------------------
    # Clustering (Leiden) - useful for grouping cells
    # -------------------------------------------------------------------------
    print(f"[{sample_name}] Computing Leiden clustering")
    sc.tl.leiden(adata, resolution=0.1)
    
    # Plot UMAP colored by clusters
    sc.pl.umap(adata, color='leiden', show=False)
    generated_plots.append(
        ("umap_leiden_clusters.png", "UMAP colored by Leiden clusters")
    )
    save_current_fig("umap_leiden_clusters.png")

    analyzed_path = None
    if save_anndata is not None:
        print(f"[{sample_name}] Saving analyzed AnnData: {save_anndata}")
        analyzed_path = write_analyzed_adata(adata, save_anndata)

    # -------------------------------------------------------------------------
    # Top velocity genes heatmap (optional but useful)
    # -------------------------------------------------------------------------
    # print(f"[{sample_name}] Computing top velocity genes")
    # scv.tl.rank_velocity_genes(adata, groupby='leiden', n_genes=50)
    # print(f"[{sample_name}] Generating top velocity genes heatmap")
    # # scv.pl.rank_velocity_genes(
    # #     adata,
    # #     n_genes=20,
    # #     sharey=False,
    # #     show=False
    # # )
    # scv.pl.scatter(adata, 
    #         basis=adata.uns["rank_velocity_genes"]["names"]["Beta"][:4])
    # print(f"[{sample_name}] Generating top velocity genes plot")
    # generated_plots.append(
    #     ("rank_velocity_genes.png",
    #      "Top velocity genes (rank_velocity_genes)")
    # )
    # save_current_fig("rank_velocity_genes.png")

    # -------------------------------------------------------------------------
    # Build simple HTML report
    # -------------------------------------------------------------------------
    print(f"[{sample_name}] Building HTML report")

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    report_filename = f"{sample_name}_scvelo_report.html"
    report_path = os.path.join(output_dir, report_filename)

    html_parts = [
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        f"<title>scVelo report - {sample_name}</title>",
        '<meta charset="utf-8" />',
        "<style>",
        "body { font-family: sans-serif; max-width: 1200px; margin: 0 auto; }",
        "h1, h2 { font-family: sans-serif; }",
        "img { max-width: 100%; height: auto; margin-bottom: 24px; }",
        ".plot-block { margin-bottom: 40px; }",
        ".subtitle { color: #555; font-size: 0.9em; }",
        "</style>",
        "</head>",
        "<body>",
        f"<h1>scVelo Report: {sample_name}</h1>",
        f'<p class="subtitle">Generated: {timestamp}</p>',
        "<h2>QC and Velocity Analysis</h2>",
    ]

    if metadata_details is not None:
        html_parts.extend(
            [
                "<h2>External cell metadata</h2>",
                "<ul>",
                f"<li>File: {escape(metadata_details['file'])}</li>",
                f"<li>Join: {escape(metadata_details['metadata_key'])} &rarr; "
                f"{escape(metadata_details['adata_key'])}</li>",
                f"<li>Matched cells: {metadata_details['matched_cells']:,}</li>",
                f"<li>Additional metadata rows ignored: "
                f"{metadata_details['extra_metadata_rows']:,}</li>",
                "</ul>",
            ]
        )

    if subset_details is not None:
        html_parts.extend(
            [
                "<h2>Cell subset</h2>",
                "<ul>",
                f"<li>Column: {escape(subset_details['column'])}</li>",
                f"<li>Values: "
                f"{escape(', '.join(subset_details['values']))}</li>",
                f"<li>Retained cells: {subset_details['retained_cells']:,} of "
                f"{subset_details['original_cells']:,}</li>",
                "</ul>",
            ]
        )

    if analyzed_path is not None:
        html_parts.extend(
            [
                "<h2>Analysis artifacts</h2>",
                "<ul>",
                f"<li>Analyzed AnnData: {escape(analyzed_path)}</li>",
                "</ul>",
            ]
        )

    for filename, title in generated_plots:
        html_parts.append('<div class="plot-block">')
        html_parts.append(f"<h3>{escape(title)}</h3>")
        html_parts.append(
            f'<img src="{escape(filename)}" alt="{escape(title)}">'
        )
        html_parts.append("</div>")

    html_parts.extend(["</body>", "</html>"])

    with open(report_path, "w") as f:
        f.write("\n".join(html_parts))

    print(f"[{sample_name}] Report written to: {report_path}")
    return report_path
