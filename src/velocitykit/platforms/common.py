"""Generic dual-run subtraction implementation for velocity-kit.

This module provides a generic implementation that works for any platform
using the dual-run subtraction method (total - exonic = unspliced).

Redesigned to use scanpy for robust data loading and AnnData operations.
"""

import logging
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

try:
    import scanpy as sc
    HAS_SCANPY = True
except ImportError:
    HAS_SCANPY = False

try:
    from tqdm.auto import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

logger = logging.getLogger(__name__)


def run_dual_subtraction(
    args,
    platform_name: str,
    total_help_text: str,
    exonic_help_text: str,
    subdirectory: str = None
):
    """
    Generic dual-run subtraction pipeline using scanpy for robust data loading.
    
    This implementation:
    1. Loads both matrices using scanpy's battle-tested read_10x_mtx()
    2. Subsets exonic matrix to match filtered cells from total matrix
    3. Aligns genes across both matrices (union of genes)
    4. Computes unspliced = total - exonic
    5. Creates AnnData with spliced/unspliced layers
    
    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments with:
        - total: Path to total counts directory
        - exonic: Path to exonic counts directory
        - genes_col: Column index for gene IDs (0=gene_ids, 1=gene_symbols)
        - out_h5ad: Output H5AD path (optional)
        - out_loom: Output loom path (optional)
    platform_name : str
        Name of the platform for logging (e.g., "PIPseq", "10x Genomics")
    total_help_text : str
        Description of the total counts input for error messages
    exonic_help_text : str
        Description of the exonic counts input for error messages
    subdirectory : str, optional
        If provided, will look for this subdirectory (e.g., "raw_feature_bc_matrix")
        when the matrix files are not found directly
        
    Raises
    ------
    ValueError
        If no output format is specified
    NotADirectoryError
        If input directories don't exist
    ImportError
        If scanpy is not installed
    """
    if not HAS_SCANPY:
        raise ImportError(
            "scanpy is required for this functionality. "
            "Install it with: pip install scanpy"
        )
    
    # Validate output arguments
    if not args.out_h5ad and not args.out_loom:
        raise ValueError(
            "At least one output format must be specified: --out-h5ad or --out-loom"
        )
    
    total_dir = Path(args.total)
    ex_dir = Path(args.exonic)

    # Directory validation
    if not total_dir.is_dir():
        raise NotADirectoryError(f"--total is not a directory: {total_dir}")
    if not ex_dir.is_dir():
        raise NotADirectoryError(f"--exonic is not a directory: {ex_dir}")

    logger.info(f"Starting {platform_name} velocity matrix construction.")
    logger.info(f"Strategy: Dual-run subtraction (unspliced = total - exonic)")
    logger.info(f"NOTE: {total_help_text}")
    logger.info(f"      {exonic_help_text}")

    # Count steps for progress bar
    num_outputs = int(bool(args.out_h5ad)) + int(bool(args.out_loom))
    steps = 4 + num_outputs
    if HAS_TQDM:
        pbar = tqdm(total=steps, desc="Pipeline", ncols=80)
    else:
        pbar = None

    def step_done():
        if pbar is not None:
            pbar.update(1)

    # Handle optional subdirectory structure (e.g., for 10x)
    if subdirectory:
        total_dir = _find_matrix_dir(total_dir, subdirectory)
        ex_dir = _find_matrix_dir(ex_dir, subdirectory)

    # Step 1: Load total matrix (filtered cells, introns included)
    logger.info("Loading total matrix (introns included) using scanpy...")
    var_names = 'gene_ids' if args.genes_col == 0 else 'gene_symbols'
    adata_total = sc.read_10x_mtx(
        total_dir,
        var_names=var_names,
        cache=False,
        gex_only=True  # Only gene expression, skip other feature types
    )
    logger.info(
        f"Total matrix: {adata_total.n_obs} cells × {adata_total.n_vars} genes "
        f"(typically FILTERED cells)"
    )
    step_done()

    # Step 2: Load exonic matrix (unfiltered cells, exons only)
    logger.info("Loading exonic matrix (RAW/UNFILTERED) using scanpy...")
    adata_exonic = sc.read_10x_mtx(
        ex_dir,
        var_names=var_names,
        cache=False,
        gex_only=True
    )
    logger.info(
        f"Exonic matrix: {adata_exonic.n_obs} cells × {adata_exonic.n_vars} genes "
        f"(should be RAW/UNFILTERED - more cells than total)"
    )
    
    # Sanity check: exonic should have >= cells as total (it's unfiltered)
    if adata_exonic.n_obs < adata_total.n_obs:
        logger.warning(
            f"⚠️  Exonic matrix has fewer cells ({adata_exonic.n_obs}) than "
            f"total matrix ({adata_total.n_obs}). "
            f"Expected exonic to be RAW/UNFILTERED with MORE cells. "
            f"Please verify your inputs!"
        )
    step_done()

    # Step 3: Align matrices and compute spliced/unspliced
    logger.info("Aligning matrices and computing spliced/unspliced counts...")
    adata = build_velocity_adata_from_anndata(adata_total, adata_exonic)
    logger.info(
        f"Final output: {adata.n_obs} cells × {adata.n_vars} genes "
        f"with 'spliced' and 'unspliced' layers"
    )
    step_done()

    # Step 4: Quality checks
    _log_quality_metrics(adata)
    step_done()

    # Write outputs
    if args.out_h5ad:
        out_h5ad = Path(args.out_h5ad)
        out_h5ad.parent.mkdir(parents=True, exist_ok=True)
        logger.info(f"Writing H5AD to {out_h5ad}")
        adata.write_h5ad(str(out_h5ad), compression='gzip')
        step_done()

    if args.out_loom:
        out_loom = Path(args.out_loom)
        out_loom.parent.mkdir(parents=True, exist_ok=True)
        logger.info(f"Writing loom to {out_loom}")
        adata.write_loom(str(out_loom))
        step_done()

    if pbar is not None:
        pbar.close()

    logger.info("✅ Successfully built velocity-compatible files!")
    if args.out_h5ad:
        logger.info(f"   H5AD: {out_h5ad}")
    if args.out_loom:
        logger.info(f"   LOOM: {out_loom}")
    logger.info("Ready for velocity analysis with scVelo or velocyto!")


def _find_matrix_dir(base_dir: Path, subdirectory: str) -> Path:
    """
    Find the directory containing matrix files.
    
    Looks for matrix.mtx.gz in base_dir first, then in base_dir/subdirectory.
    Common for 10x data which may have filtered_feature_bc_matrix/ or
    raw_feature_bc_matrix/ subdirectories.
    
    Parameters
    ----------
    base_dir : Path
        Base directory to search
    subdirectory : str
        Subdirectory name to check if files not found in base_dir
        
    Returns
    -------
    Path
        Directory containing the matrix files
        
    Raises
    ------
    FileNotFoundError
        If matrix.mtx.gz not found in either location
    """
    # Check base directory first
    matrix_file = base_dir / "matrix.mtx.gz"
    if matrix_file.exists():
        return base_dir
    
    # Also check for uncompressed version
    matrix_file_uncompressed = base_dir / "matrix.mtx"
    if matrix_file_uncompressed.exists():
        return base_dir
    
    # Try subdirectory
    subdir_path = base_dir / subdirectory
    subdir_matrix = subdir_path / "matrix.mtx.gz"
    if subdir_matrix.exists():
        logger.info(f"Found matrix files in subdirectory: {subdir_path}")
        return subdir_path
    
    subdir_matrix_uncompressed = subdir_path / "matrix.mtx"
    if subdir_matrix_uncompressed.exists():
        logger.info(f"Found matrix files in subdirectory: {subdir_path}")
        return subdir_path
    
    raise FileNotFoundError(
        f"Could not find matrix.mtx.gz or matrix.mtx in:\n"
        f"  - {base_dir}\n"
        f"  - {subdir_path}\n"
        f"Please verify the directory structure."
    )


def build_velocity_adata_from_anndata(
    adata_total: ad.AnnData,
    adata_exonic: ad.AnnData
) -> ad.AnnData:
    """
    Build velocity AnnData by aligning and subtracting exonic from total counts.
    
    Strategy:
    1. Use cell barcodes from total matrix (filtered cells)
    2. Subset exonic matrix to only include these filtered cells
    3. Take union of genes from both matrices
    4. Align both matrices to this gene union
    5. Compute: unspliced = total - exonic, spliced = exonic
    
    Parameters
    ----------
    adata_total : ad.AnnData
        Total counts (exons + introns), typically FILTERED cells
    adata_exonic : ad.AnnData
        Exonic counts only, typically UNFILTERED cells (raw matrix)
        
    Returns
    -------
    ad.AnnData
        AnnData with:
        - X: spliced counts (exonic)
        - layers['spliced']: same as X
        - layers['unspliced']: total - exonic (intronic)
        - obs: from total matrix
        - var: union of genes
    """
    # Step 1: Get filtered cell barcodes from total matrix
    barcodes_filtered = adata_total.obs_names
    n_filtered = len(barcodes_filtered)
    logger.info(f"Using {n_filtered} filtered cells from total matrix as reference")
    
    # Step 2: Find overlap with exonic matrix
    common_cells = barcodes_filtered.intersection(adata_exonic.obs_names)
    n_common = len(common_cells)
    overlap_pct = 100 * n_common / n_filtered
    
    logger.info(
        f"Found {n_common}/{n_filtered} cells in exonic matrix "
        f"({overlap_pct:.1f}% overlap)"
    )
    
    # Quality check: should have high overlap
    if overlap_pct < 50:
        logger.error(
            f"❌ LOW OVERLAP WARNING: Only {overlap_pct:.1f}% of filtered cells "
            f"found in exonic matrix. This suggests a barcode mismatch problem!"
        )
        raise ValueError(
            "Insufficient barcode overlap between total and exonic matrices. "
            "Please verify that both matrices are from the same sample and "
            "use compatible barcode formats."
        )
    elif overlap_pct < 90:
        logger.warning(
            f"⚠️  Moderate overlap: {overlap_pct:.1f}%. "
            f"Missing {n_filtered - n_common} cells from exonic matrix."
        )
    
    # Step 3: Subset both matrices to common cells only
    logger.info("Subsetting matrices to common cells...")
    adata_total = adata_total[common_cells, :].copy()
    adata_exonic = adata_exonic[common_cells, :].copy()
    
    # Step 4: Gene alignment - take union of genes
    genes_total = set(adata_total.var_names)
    genes_exonic = set(adata_exonic.var_names)
    genes_union = sorted(genes_total | genes_exonic)
    genes_common = genes_total & genes_exonic
    
    logger.info(
        f"Gene statistics:\n"
        f"  Total matrix:  {len(genes_total)} genes\n"
        f"  Exonic matrix: {len(genes_exonic)} genes\n"
        f"  Common genes:  {len(genes_common)} genes\n"
        f"  Union (final): {len(genes_union)} genes"
    )
    
    # Step 5: Reindex both matrices to gene union
    logger.info("Aligning genes to union...")
    adata_total_aligned = _reindex_genes(adata_total, genes_union)
    adata_exonic_aligned = _reindex_genes(adata_exonic, genes_union)
    
    # Step 6: Compute unspliced = total - exonic
    logger.info("Computing unspliced counts (total - exonic)...")
    X_unspliced = adata_total_aligned.X - adata_exonic_aligned.X
    
    # Handle negative values (can occur due to sampling/technical noise)
    n_negative = _clip_negative_values(X_unspliced)
    if n_negative > 0:
        total_entries = X_unspliced.shape[0] * X_unspliced.shape[1]
        pct_negative = 100 * n_negative / total_entries
        logger.warning(
            f"Clipped {n_negative:,} negative values to 0 "
            f"({pct_negative:.3f}% of entries)"
        )
        if pct_negative > 5:
            logger.warning(
                f"⚠️  HIGH NEGATIVE RATE: {pct_negative:.1f}% negative values! "
                f"This may indicate quality issues with the input data."
            )
    
    # Step 7: Build final AnnData with spliced/unspliced layers
    logger.info("Building final AnnData object...")
    adata = ad.AnnData(
        X=adata_exonic_aligned.X.copy(),  # Main matrix = spliced (exonic)
        obs=adata_total.obs.copy(),
        var=pd.DataFrame(index=genes_union)
    )
    
    adata.layers['spliced'] = adata_exonic_aligned.X.copy()
    adata.layers['unspliced'] = X_unspliced
    
    # Format obs/var to match velocyto standard:
    # Column attributes (cells): CellID
    # Row attributes (genes): Gene, Accession
    
    # Ensure indices are strings
    adata.obs.index = adata.obs.index.astype(str)
    adata.var.index = adata.var.index.astype(str)
    
    # Add CellID column (velocyto standard for cell barcodes)
    if 'CellID' not in adata.obs.columns:
        adata.obs['CellID'] = adata.obs.index.astype(str)
    
    # Add Gene column (velocyto standard for gene names/symbols)
    # Assumes var_names are gene symbols or IDs
    if 'Gene' not in adata.var.columns:
        adata.var['Gene'] = adata.var.index.astype(str)
    
    # Add Accession column (velocyto standard for gene IDs like ENSMUSG...)
    # If var_names look like Ensembl IDs, use them; otherwise duplicate Gene
    if 'Accession' not in adata.var.columns:
        adata.var['Accession'] = adata.var.index.astype(str)
    
    # Make sure names are unique to avoid write issues
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    
    # Add useful metadata
    adata.uns['velocity_params'] = {
        'method': 'dual_run_subtraction',
        'n_cells_total': n_filtered,
        'n_cells_overlap': n_common,
        'overlap_pct': overlap_pct,
        'n_genes_total': len(genes_total),
        'n_genes_exonic': len(genes_exonic),
        'n_genes_union': len(genes_union),
        'n_negative_clipped': int(n_negative)
    }
    
    return adata


def _reindex_genes(adata: ad.AnnData, genes_union: list) -> ad.AnnData:
    """
    Reindex AnnData to have all genes in genes_union.
    
    Missing genes are filled with zeros. This is memory-efficient using
    sparse matrices.
    
    Parameters
    ----------
    adata : ad.AnnData
        Input AnnData
    genes_union : list
        List of all genes to include
        
    Returns
    -------
    ad.AnnData
        Reindexed AnnData with all genes from genes_union
    """
    current_genes = list(adata.var_names)
    missing_genes = [g for g in genes_union if g not in current_genes]
    
    if not missing_genes:
        # Already has all genes, just reorder
        return adata[:, genes_union].copy()
    
    logger.debug(f"Adding {len(missing_genes)} missing genes (filled with zeros)")
    
    # Create zero matrix for missing genes (use sparse for efficiency)
    n_obs = adata.n_obs
    n_missing = len(missing_genes)
    
    if sp.issparse(adata.X):
        zero_matrix = sp.csr_matrix((n_obs, n_missing), dtype=adata.X.dtype)
    else:
        zero_matrix = np.zeros((n_obs, n_missing), dtype=adata.X.dtype)
    
    # Create AnnData for missing genes
    missing_var = pd.DataFrame(index=missing_genes)
    adata_missing = ad.AnnData(
        X=zero_matrix,
        obs=adata.obs,
        var=missing_var
    )
    
    # Concatenate along gene axis (axis=1)
    adata_concat = ad.concat([adata, adata_missing], axis=1, join='outer')
    
    # Reorder to match genes_union
    return adata_concat[:, genes_union].copy()


def _clip_negative_values(X) -> int:
    """
    Clip negative values in matrix to 0 (in-place).
    
    Parameters
    ----------
    X : array-like or sparse matrix
        Input matrix
        
    Returns
    -------
    int
        Number of negative values clipped
    """
    if sp.issparse(X):
        # Sparse matrix
        negative_mask = X.data < 0
        n_negative = negative_mask.sum()
        if n_negative > 0:
            X.data[negative_mask] = 0
            X.eliminate_zeros()
    else:
        # Dense matrix
        negative_mask = X < 0
        n_negative = negative_mask.sum()
        if n_negative > 0:
            X[negative_mask] = 0
    
    return int(n_negative)


def _log_quality_metrics(adata: ad.AnnData):
    """Log quality metrics for the velocity data."""
    X_spliced = adata.layers['spliced']
    X_unspliced = adata.layers['unspliced']
    
    # Compute sparsity
    if sp.issparse(X_spliced):
        spliced_nonzero = X_spliced.nnz
        unspliced_nonzero = X_unspliced.nnz
    else:
        spliced_nonzero = np.count_nonzero(X_spliced)
        unspliced_nonzero = np.count_nonzero(X_unspliced)
    
    total_entries = adata.n_obs * adata.n_vars
    spliced_density = 100 * spliced_nonzero / total_entries
    unspliced_density = 100 * unspliced_nonzero / total_entries
    
    # Compute count statistics
    if sp.issparse(X_spliced):
        spliced_total = X_spliced.data.sum()
        unspliced_total = X_unspliced.data.sum()
    else:
        spliced_total = X_spliced.sum()
        unspliced_total = X_unspliced.sum()
    
    logger.info(
        f"\n📊 Quality Metrics:\n"
        f"  Spliced counts:   {spliced_total:,.0f} ({spliced_density:.2f}% density)\n"
        f"  Unspliced counts: {unspliced_total:,.0f} ({unspliced_density:.2f}% density)\n"
        f"  Unspliced ratio:  {100 * unspliced_total / (spliced_total + unspliced_total):.1f}%"
    )


def add_standard_arguments(parser, platform_name: str, default_genes_col: int = 0):
    """
    Add standard arguments for dual-run subtraction commands.
    
    Parameters
    ----------
    parser : argparse.ArgumentParser
        Argument parser to add arguments to
    platform_name : str
        Platform name for help text (e.g., "PIPseq", "10x Genomics")
    default_genes_col : int
        Default column index for gene IDs in features.tsv
        0 = gene IDs, 1 = gene symbols
    """
    parser.add_argument(
        "--total",
        required=True,
        help=(
            f"Directory with {platform_name} total counts (introns + exons). "
            f"Should contain matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz. "
            f"Typically uses FILTERED cells."
        ),
    )
    parser.add_argument(
        "--exonic",
        required=True,
        help=(
            f"Directory with {platform_name} exonic counts (exons only). "
            f"Should contain matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz. "
            f"MUST use RAW/UNFILTERED count matrix (before cell calling)."
        ),
    )
    parser.add_argument(
        "--genes-col",
        type=int,
        default=default_genes_col,
        help=(
            f"Column index in features.tsv to use as gene identifier. "
            f"0 = gene IDs (e.g., ENSMUSG...), 1 = gene symbols (e.g., Actb). "
            f"Default: {default_genes_col}"
        ),
    )
    parser.add_argument(
        "--out-h5ad",
        help="Output .h5ad file path (optional if --out-loom is provided).",
    )
    parser.add_argument(
        "--out-loom",
        help="Output .loom file path (optional if --out-h5ad is provided).",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=1,
        help="Increase verbosity level. Use -v for INFO, -vv for DEBUG.",
    )
