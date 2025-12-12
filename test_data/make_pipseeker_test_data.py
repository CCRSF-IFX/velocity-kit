#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def write_gzip_tsv(path, df):
    """Write a small TSV file compressed with gzip."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as f:
        df.to_csv(f, sep="\t", header=False, index=False)


def write_10x_mtx_dir(out_dir, X, barcodes, genes):
    """
    Write a 10x-like matrix directory:
      - matrix.mtx.gz
      - barcodes.tsv.gz
      - features.tsv.gz
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Matrix (genes x barcodes)
    mtx_path = out_dir / "matrix.mtx.gz"
    with gzip.open(mtx_path, "wb") as f:
        scipy.io.mmwrite(f, X)

    # Barcodes
    barcodes_path = out_dir / "barcodes.tsv.gz"
    write_gzip_tsv(barcodes_path, barcodes.to_frame())

    # Features: simple 3-column features.tsv: gene_id, gene_name, feature_type
    features_df = pd.DataFrame(
        {
            0: genes,                 # gene_id
            1: genes,                 # gene_name (same for test)
            2: "Gene Expression",     # feature type (10x style)
        }
    )
    features_path = out_dir / "features.tsv.gz"
    write_gzip_tsv(features_path, features_df)


def generate_synthetic_counts(
    n_genes,
    n_cells_filtered,
    n_cells_raw,
    lambda_total,
    lambda_exonic_fraction,
    seed,
):
    """
    Generate synthetic total and exonic count matrices with:

        total = exonic + intronic (FILTERED cells)
        exonic = exonic-only (RAW/UNFILTERED cells, includes low-quality cells)

    The exonic matrix has MORE barcodes than total (mimics real PIPseq workflow).
    For the overlapping cells:
        exonic <= total elementwise

    total ~ Poisson(lambda_total) per gene-cell
    exonic ~ Binomial(total, p_exon) for filtered cells
    exonic ~ Poisson(lambda_total * p_exon * 0.1) for extra raw cells (background)
    p_exon varies by gene around lambda_exonic_fraction.
    """
    rng = np.random.default_rng(seed)

    # Gene IDs
    genes = pd.Series(["Gene%05d" % i for i in range(n_genes)], dtype=str)
    
    # Barcodes: filtered cells are a subset of raw cells
    barcodes_raw = pd.Series(["BC%05d" % i for i in range(n_cells_raw)], dtype=str)
    barcodes_filtered = barcodes_raw[:n_cells_filtered].copy()

    # Total counts (genes x filtered_cells) - this is the FILTERED dataset
    total_counts = rng.poisson(lam=lambda_total, size=(n_genes, n_cells_filtered))

    # Gene-specific exon fractions, with some variation
    exon_frac_per_gene = np.clip(
        rng.normal(loc=lambda_exonic_fraction, scale=0.1, size=n_genes),
        0.05,
        0.95,
    )
    
    # Exonic counts for FILTERED cells: binomial(total, p_exon)
    exon_frac_matrix_filtered = exon_frac_per_gene[:, None]
    exonic_counts_filtered = rng.binomial(n=total_counts, p=exon_frac_matrix_filtered)
    assert np.all(exonic_counts_filtered <= total_counts)

    # Exonic counts for EXTRA RAW cells (background/low-quality cells)
    # These have lower counts (10% of typical exonic counts)
    n_extra_cells = n_cells_raw - n_cells_filtered
    lambda_background = lambda_total * lambda_exonic_fraction * 0.1
    exonic_counts_background = rng.poisson(
        lam=lambda_background, 
        size=(n_genes, n_extra_cells)
    )
    
    # Combine: exonic_raw = [filtered cells | background cells]
    exonic_counts_raw = np.hstack([exonic_counts_filtered, exonic_counts_background])

    X_total = sp.coo_matrix(total_counts)
    X_exonic = sp.coo_matrix(exonic_counts_raw)

    return X_total, X_exonic, genes, barcodes_filtered, barcodes_raw


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate synthetic 10x-like PIPseeker test data:\n"
            "  - total/       (total = exonic + intronic)\n"
            "  - exonic_raw/  (exonic-only, mimics RAW/UNFILTERED exons-only run)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--out-root",
        required=True,
        help="Root output directory for test data.",
    )
    parser.add_argument(
        "--n-genes",
        type=int,
        default=1000,
        help="Number of genes (default: 1000).",
    )
    parser.add_argument(
        "--n-cells-filtered",
        type=int,
        default=200,
        help="Number of FILTERED cells/barcodes in total run (default: 200).",
    )
    parser.add_argument(
        "--n-cells-raw",
        type=int,
        default=500,
        help="Number of RAW/UNFILTERED cells/barcodes in exonic run (default: 500, must be >= n-cells-filtered).",
    )
    parser.add_argument(
        "--lambda-total",
        type=float,
        default=1.0,
        help="Mean total counts per gene-cell (Poisson rate; default: 1.0).",
    )
    parser.add_argument(
        "--lambda-exonic-fraction",
        type=float,
        default=0.7,
        help="Average fraction of total counts that are exonic (default: 0.7).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Random seed for reproducibility (default: 1).",
    )

    args = parser.parse_args()

    # Validate arguments
    if args.n_cells_raw < args.n_cells_filtered:
        parser.error(
            f"--n-cells-raw ({args.n_cells_raw}) must be >= --n-cells-filtered ({args.n_cells_filtered})"
        )

    out_root = Path(args.out_root)
    out_total = out_root / "total"
    out_exonic = out_root / "exonic_raw"

    print(
        "Generating synthetic PIPseeker-style test data with parameters:\n"
        "  genes                    = {n_genes}\n"
        "  cells (filtered, total)  = {n_cells_filtered}\n"
        "  cells (raw, exonic)      = {n_cells_raw}\n"
        "  lambda_total             = {lambda_total}\n"
        "  lambda_exonic_fraction   = {lambda_exonic_fraction}\n"
        "  seed                     = {seed}\n"
        "  out_root                 = {out_root}\n".format(
            n_genes=args.n_genes,
            n_cells_filtered=args.n_cells_filtered,
            n_cells_raw=args.n_cells_raw,
            lambda_total=args.lambda_total,
            lambda_exonic_fraction=args.lambda_exonic_fraction,
            seed=args.seed,
            out_root=out_root,
        )
    )

    X_total, X_exonic, genes, barcodes_filtered, barcodes_raw = generate_synthetic_counts(
        n_genes=args.n_genes,
        n_cells_filtered=args.n_cells_filtered,
        n_cells_raw=args.n_cells_raw,
        lambda_total=args.lambda_total,
        lambda_exonic_fraction=args.lambda_exonic_fraction,
        seed=args.seed,
    )

    print("Writing total (exonic + intronic) FILTERED dataset → total/ ...")
    write_10x_mtx_dir(out_total, X_total, barcodes_filtered, genes)

    print("Writing exonic-only RAW/UNFILTERED dataset → exonic_raw/ ...")
    write_10x_mtx_dir(out_exonic, X_exonic, barcodes_raw, genes)

    print("Done.")
    print("  Total dir:  %s" % out_total)
    print("  Exonic dir: %s" % out_exonic)


if __name__ == "__main__":
    main()
