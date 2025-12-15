"""Core functionality for pipseq-velocity package."""

import gzip
import logging
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad

logger = logging.getLogger(__name__)


def _open_maybe_gzip(path: str):
    """Open a file, handling gzipped files automatically."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def load_10x_mtx(
    matrix_path: Path,
    barcodes_path: Path,
    features_path: Path,
    genes_col: int = 0
) -> Tuple[sp.csr_matrix, pd.Series, pd.Series]:
    """
    Load a 10x-like MTX matrix (genes x barcodes) + barcodes + features.

    Parameters
    ----------
    matrix_path : Path
        Path to the matrix.mtx(.gz) file
    barcodes_path : Path
        Path to the barcodes.tsv(.gz) file
    features_path : Path
        Path to the features.tsv(.gz) file
    genes_col : int, optional
        Column index in features.tsv to use as gene ID (default: 0)

    Returns
    -------
    X : scipy.sparse.csr_matrix
        The count matrix (genes x barcodes)
    barcodes : pd.Series
        The barcode IDs
    genes : pd.Series
        The gene IDs

    Raises
    ------
    FileNotFoundError
        If any required file does not exist
    ValueError
        If dimensions don't match or genes_col is out of range
    """
    matrix_path = Path(matrix_path)
    barcodes_path = Path(barcodes_path)
    features_path = Path(features_path)

    # Argument validation: existence
    for p in [matrix_path, barcodes_path, features_path]:
        if not p.exists():
            raise FileNotFoundError(f"Required file not found: {p}")

    logger.info(f"Reading matrix: {matrix_path}")
    X = scipy.io.mmread(str(matrix_path)).tocsr()
    if X.nnz == 0:
        logger.warning(f"Matrix {matrix_path} has 0 non-zero entries.")

    logger.info(f"Reading barcodes: {barcodes_path}")
    with _open_maybe_gzip(barcodes_path) as f:
        barcodes = pd.read_csv(f, header=None, sep="\t")[0].astype(str)

    logger.info(f"Reading features: {features_path}")
    with _open_maybe_gzip(features_path) as f:
        features = pd.read_csv(f, header=None, sep="\t")

    if genes_col < 0 or genes_col >= features.shape[1]:
        raise ValueError(
            f"genes_col {genes_col} is out of range for features.tsv "
            f"(n_cols = {features.shape[1]})"
        )

    genes = features[genes_col].astype(str)

    # Dimension checks
    if X.shape[0] != len(genes):
        raise ValueError(
            f"Matrix gene dimension {X.shape[0]} != number of genes {len(genes)} "
            f"in {features_path}"
        )
    if X.shape[1] != len(barcodes):
        raise ValueError(
            f"Matrix cell dimension {X.shape[1]} != number of barcodes {len(barcodes)} "
            f"in {barcodes_path}"
        )

    # Simple format checks
    if genes.duplicated().any():
        logger.warning("Gene IDs contain duplicates.")
    if barcodes.duplicated().any():
        logger.warning("Barcodes contain duplicates.")
    if (genes.str.len() == 0).any():
        logger.warning("Some gene IDs are empty strings.")
    if (barcodes.str.len() == 0).any():
        logger.warning("Some barcodes are empty strings.")

    logger.info(
        f"Loaded matrix with shape {X.shape[0]} genes x {X.shape[1]} barcodes."
    )

    return X, barcodes, genes


def align_and_union(
    X1: sp.csr_matrix,
    bc1: pd.Series,
    g1: pd.Series,
    X2: sp.csr_matrix,
    bc2: pd.Series,
    g2: pd.Series
) -> Tuple[sp.csr_matrix, sp.csr_matrix, pd.Index, pd.Index]:
    """
    Align two matrices to the union of genes but only the barcodes from X1 (total/filtered).
    
    This matches the typical velocity workflow where:
    - X1 (total) contains FILTERED cells
    - X2 (exonic) contains RAW/UNFILTERED cells (superset of X1 barcodes)
    - Output matrices contain only the filtered cells from X1
    
    Parameters
    ----------
    X1 : scipy.sparse.csr_matrix
        Total matrix (genes x barcodes) - typically FILTERED cells
    bc1 : pd.Series
        Barcodes for total matrix (filtered cells to keep)
    g1 : pd.Series
        Genes for total matrix
    X2 : scipy.sparse.csr_matrix
        Exonic matrix (genes x barcodes) - typically RAW/UNFILTERED cells
    bc2 : pd.Series
        Barcodes for exonic matrix (should be superset of bc1)
    g2 : pd.Series
        Genes for exonic matrix

    Returns
    -------
    X1u : scipy.sparse.csr_matrix
        Total matrix aligned to union genes and bc1 barcodes
    X2u : scipy.sparse.csr_matrix
        Exonic matrix aligned to union genes and bc1 barcodes
    genes_union : pd.Index
        Union of gene IDs
    barcodes_final : pd.Index
        Barcodes from X1 (filtered cells only)
    """
    logger.info("Aligning matrices...")
    
    # Use union of genes (keep all genes from both runs)
    genes_union = pd.Index(sorted(set(g1) | set(g2)))
    
    # Use barcodes from X1 only (filtered cells)
    barcodes_final = pd.Index(sorted(bc1))
    
    logger.info(
        f"Output dimensions: {len(genes_union)} genes x {len(barcodes_final)} barcodes (from total/filtered)"
    )
    
    # Check barcode overlap
    set_bc1 = set(bc1)
    set_bc2 = set(bc2)
    barcodes_intersection = set_bc1 & set_bc2
    
    logger.info(
        f"Barcode overlap: {len(barcodes_intersection)}/{len(bc1)} "
        f"({100*len(barcodes_intersection)/len(bc1):.1f}%) of filtered cells found in exonic run"
    )
    
    if len(barcodes_intersection) < len(bc1) * 0.5:
        logger.warning(
            f"Only {len(barcodes_intersection)}/{len(bc1)} filtered barcodes found in exonic run. "
            f"This may indicate a barcode mismatch between runs."
        )
    
    # Check if exonic has MORE barcodes (expected for raw/unfiltered)
    if len(bc2) > len(bc1):
        logger.info(
            f"Exonic run has {len(bc2)} barcodes (raw/unfiltered), "
            f"total run has {len(bc1)} barcodes (filtered). This is expected."
        )
    elif len(bc2) == len(bc1):
        logger.warning(
            "Exonic and total runs have the same number of barcodes. "
            "For velocity analysis, the exonic run should typically be RAW/UNFILTERED "
            "(more barcodes than the filtered total run)."
        )
    else:
        logger.warning(
            "Exonic run has FEWER barcodes than total run. "
            "This is unusual - make sure you're using the RAW/UNFILTERED exonic matrix."
        )
    
    # Create index maps
    g1_map = pd.Series(np.arange(len(g1)), index=g1)
    g2_map = pd.Series(np.arange(len(g2)), index=g2)
    bc1_map = pd.Series(np.arange(len(bc1)), index=bc1)
    bc2_map = pd.Series(np.arange(len(bc2)), index=bc2)

    def scatter(X_src, gmap, bcmap, label):
        """Align source matrix to final gene and barcode space."""
        genes_common = genes_union.intersection(gmap.index)
        bc_common = barcodes_final.intersection(bcmap.index)

        logger.info(
            f"[{label}] Overlap: {len(genes_common)} genes, {len(bc_common)} barcodes."
        )

        if len(genes_common) == 0 or len(bc_common) == 0:
            logger.warning(
                f"[{label}] No overlapping genes or barcodes. Matrix will be all zeros."
            )
            return sp.csr_matrix((len(genes_union), len(barcodes_final)))

        # Get source indices
        src_g_idx = gmap[genes_common].values
        src_bc_idx = bcmap[bc_common].values
        X_sub = X_src[src_g_idx[:, None], src_bc_idx]

        # Get destination indices
        dest_g_idx = genes_union.get_indexer(genes_common)
        dest_bc_idx = barcodes_final.get_indexer(bc_common)

        # Remap to destination space
        X_sub = X_sub.tocoo()
        X_sub.row = dest_g_idx[X_sub.row]
        X_sub.col = dest_bc_idx[X_sub.col]

        X_aligned = sp.csr_matrix(
            (X_sub.data, (X_sub.row, X_sub.col)),
            shape=(len(genes_union), len(barcodes_final)),
        )
        return X_aligned

    X1u = scatter(X1, g1_map, bc1_map, label="total")
    X2u = scatter(X2, g2_map, bc2_map, label="exonic")

    return X1u, X2u, genes_union, barcodes_final


def build_velocity_adata(
    X_total: sp.csr_matrix,
    X_exon: sp.csr_matrix,
    genes: pd.Index,
    barcodes: pd.Index
) -> ad.AnnData:
    """
    Build an AnnData object with spliced/unspliced layers for velocity analysis.

    Parameters
    ----------
    X_total : scipy.sparse.csr_matrix
        Total counts matrix (genes x barcodes)
    X_exon : scipy.sparse.csr_matrix
        Exonic counts matrix (genes x barcodes)
    genes : pd.Index
        Gene IDs
    barcodes : pd.Index
        Barcode IDs

    Returns
    -------
    adata : anndata.AnnData
        AnnData object with spliced/unspliced layers
    """
    logger.info("Building spliced/unspliced layers...")
    X_spliced = X_exon.copy().tocsr()
    X_unspliced = (X_total - X_spliced).tocsr()

    # Clip negatives
    if X_unspliced.nnz > 0:
        neg_mask = X_unspliced.data < 0
        if np.any(neg_mask):
            logger.warning(
                f"Found {neg_mask.sum()} negative entries in unspliced matrix; "
                f"clipping to 0."
            )
            X_unspliced.data[neg_mask] = 0
            X_unspliced.eliminate_zeros()

    logger.info("Creating AnnData object (cells x genes)...")
    adata = ad.AnnData(
        X=X_spliced.T,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=genes),
    )

    adata.layers["spliced"] = X_spliced.T
    adata.layers["unspliced"] = X_unspliced.T

    return adata


def run_scvelo_preprocessing(adata: ad.AnnData) -> None:
    """
    Run basic scVelo preprocessing to prepare for velocity analysis.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with spliced/unspliced layers

    Raises
    ------
    ImportError
        If scvelo is not installed
    """
    logger.info("Running scVelo preprocessing...")

    try:
        import scvelo as scv
    except ImportError:
        logger.error(
            "scvelo is not installed but preprocessing was requested. "
            "Install scvelo or run without this option."
        )
        raise

    scv.logging.print_version_and_date()
    # scv.pp.filter_genes(adata, min_shared_counts=10)
    # scv.pp.normalize_per_cell(adata)
    # scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
    # scv.pp.log1p(adata)
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.pl.proportions(adata)
    
    logger.info("scVelo preprocessing finished.")
