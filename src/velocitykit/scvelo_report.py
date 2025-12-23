"""Generate scVelo analysis reports from loom files."""

from __future__ import annotations

import os
from datetime import datetime
from typing import Optional

import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt

def run_scvelo_and_generate_report(
    loom_path: str,
    output_dir: str,
    sample_name: Optional[str] = None,
) -> str:
    """
    Run a standard scVelo pipeline on a loom file and generate an HTML report
    with QC and analysis plots.

    Parameters
    ----------
    loom_path : str
        Path to the input .loom file.
    output_dir : str
        Directory where plots and HTML report will be written.
    sample_name : str, optional
        Name used in plot titles and report filename. If None, derived from loom_path.

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
        sample_name = os.path.splitext(os.path.basename(loom_path))[0]

    # make scvelo write figures into output_dir
    scv.settings.figdir = output_dir
    scv.settings.set_figure_params(frameon=False, dpi=120)

    # -------------------------------------------------------------------------
    # Load data
    # -------------------------------------------------------------------------
    print(f"[{sample_name}] Reading loom file: {loom_path}")
    adata = sc.read(loom_path)
    
    # Fix categorical columns that may cause issues with newer pandas
    # Convert any categorical columns to regular strings to avoid
    # "property 'categories' of 'Categorical' object has no setter" errors
    for col in adata.obs.columns:
        if hasattr(adata.obs[col], 'cat'):
            adata.obs[col] = adata.obs[col].astype(str)
    for col in adata.var.columns:
        if hasattr(adata.var[col], 'cat'):
            adata.var[col] = adata.var[col].astype(str)

    print(f"[{sample_name}] Data loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # -------------------------------------------------------------------------
    # scVelo preprocessing - following scvelo_example.py
    # -------------------------------------------------------------------------
    print(f"[{sample_name}] Running preprocessing and filtering")
    
    # Filter and normalize (combined step from scVelo)
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    
    # Compute moments for velocity estimation
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

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
    
    # Compute PCA and neighbors
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
    sc.tl.umap(adata)

    # -------------------------------------------------------------------------
    # Optional: Additional metrics for QC
    # -------------------------------------------------------------------------
    print(f"[{sample_name}] Computing QC metrics")
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    
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

    # 2. Number of genes per cell
    sc.pl.violin(
        adata,
        ["n_genes_by_counts"],
        jitter=0.4,
        multi_panel=False,
        show=False
    )
    generated_plots.append(
        ("qc_n_genes_by_counts.png", "QC: Number of genes per cell")
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

    for filename, title in generated_plots:
        html_parts.append('<div class="plot-block">')
        html_parts.append(f"<h3>{title}</h3>")
        html_parts.append(f'<img src="{filename}" alt="{title}">')
        html_parts.append("</div>")

    html_parts.extend(["</body>", "</html>"])

    with open(report_path, "w") as f:
        f.write("\n".join(html_parts))

    print(f"[{sample_name}] Report written to: {report_path}")
    return report_path
