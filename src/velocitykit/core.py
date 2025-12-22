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


from typing import Tuple
import numpy as np
import pandas as pd
import scipy.sparse as sp
import logging

logger = logging.getLogger(__name__)

def align_and_union(
    X1: sp.csr_matrix,
    bc1: pd.Series,
    g1: pd.Series,
    X2: sp.csr_matrix,
    bc2: pd.Series,
    g2: pd.Series,
) -> Tuple[sp.csr_matrix, sp.csr_matrix, pd.Index, pd.Index]:
    """
    Align two matrices to the SAME gene list/order (must match exactly) and to the
    barcodes from X1 (filtered cells).

    Differences from the previous strict version:
      - Duplicate gene labels are ALLOWED and PRESERVED, as long as g1 and g2
        are exactly equal as sequences (same labels, same order).
      - Barcodes must remain unique.
    """

    # Validate shapes
    if X1.shape[0] != len(g1):
        raise ValueError(f"X1 has {X1.shape[0]} rows but g1 has length {len(g1)}")
    if X1.shape[1] != len(bc1):
        raise ValueError(f"X1 has {X1.shape[1]} columns but bc1 has length {len(bc1)}")
    if X2.shape[0] != len(g2):
        raise ValueError(f"X2 has {X2.shape[0]} rows but g2 has length {len(g2)}")
    if X2.shape[1] != len(bc2):
        raise ValueError(f"X2 has {X2.shape[1]} columns but bc2 has length {len(bc2)}")

    # Barcode uniqueness enforced
    if bc1.duplicated().any():
        raise ValueError("Duplicate barcode labels found in bc1 (total). Aborting.")
    if bc2.duplicated().any():
        raise ValueError("Duplicate barcode labels found in bc2 (exonic). Aborting.")

    # Genes: allow duplicates, but require exact sequence equality
    if len(g1) != len(g2) or not all(a == b for a, b in zip(g1, g2)):
        # helpful error showing first few differences
        n_show = 8
        preview = []
        for i, (a, b) in enumerate(zip(g1, g2)):
            if a != b:
                preview.append((i, a, b))
            if len(preview) >= n_show:
                break
        logger.error(
            "Gene lists do not match exactly between total and exonic matrices. "
            "When duplicates are allowed, the gene *sequences* must be identical "
            "(same labels in the same order). Aborting."
        )
        if preview:
            logger.error("First mismatches (index, total_gene, exonic_gene): %s", preview)
        raise ValueError("Gene lists and/or order do not match between X1 and X2. Aborting.")

    # Enforce bc1 ⊆ bc2
    set_bc1 = set(bc1)
    set_bc2 = set(bc2)
    if not set_bc1.issubset(set_bc2):
        missing = sorted(list(set_bc1 - set_bc2))[:20]
        logger.error(
            "The following filtered barcodes (from bc1) were NOT found in the exonic run bc2 "
            "(showing up to 20 missing): %s", missing
        )
        raise ValueError("Some bc1 barcodes are not present in bc2. Aborting.")

    # Final outputs: genes_final uses g1 order, barcodes_final uses bc1 order
    genes_final = pd.Index(g1)
    barcodes_final = pd.Index(bc1)

    logger.info("Building aligned matrices with %d genes (order preserved) and %d barcodes (from total).",
                len(genes_final), len(barcodes_final))

    # For genes: because g1 == g2 position-wise, the row mapping from X2 -> output is identity:
    # src_row_i in X2 should map to dest_row_i in output (same i)
    n_genes = len(genes_final)
    # But to be robust, build explicit row_map arrays (global -> dest)
    row_map_X2 = np.arange(n_genes, dtype=int)   # X2 row i -> dest row i
    # For X1 same identity mapping
    row_map_X1 = np.arange(n_genes, dtype=int)

    # For columns: map bc2 (superset) to dest indices given by bc1
    # Build dict: barcode -> dest_col_index
    barcode_to_dest_idx = {b: i for i, b in enumerate(barcodes_final)}

    # Build col_map for X2 (size = X2.shape[1]); cols not in bc1 get -1 (we drop them)
    col_map_X2 = np.full(X2.shape[1], -1, dtype=int)
    # Map each bc2 position to dest index if present in bc1
    bc2_to_src_idx = {b: i for i, b in enumerate(bc2)}
    for b, src_j in bc2_to_src_idx.items():
        if b in barcode_to_dest_idx:
            col_map_X2[src_j] = barcode_to_dest_idx[b]
        # else stays -1 (we drop raw-only barcodes)

    # Build col_map for X1 (should map all cols because bc1 defines dest columns)
    bc1_to_src_idx = {b: i for i, b in enumerate(bc1)}
    col_map_X1 = np.full(X1.shape[1], -1, dtype=int)
    for b, src_j in bc1_to_src_idx.items():
        col_map_X1[src_j] = barcode_to_dest_idx[b]

    # Remap X2 entries (COO -> filter -> remap)
    X2_coo = X2.tocoo()
    r2 = X2_coo.row
    c2 = X2_coo.col
    # Keep only entries where column maps into bc1 (col_map_X2 != -1)
    keep2 = (col_map_X2[c2] != -1)
    if not np.any(keep2):
        logger.warning("No overlapping entries found when aligning X2; returning zero matrix.")
        X2u = sp.csr_matrix((n_genes, len(barcodes_final)))
    else:
        new_row2 = row_map_X2[r2[keep2]]   # positional mapping of rows
        new_col2 = col_map_X2[c2[keep2]]
        new_data2 = X2_coo.data[keep2]
        X2u = sp.csr_matrix((new_data2, (new_row2, new_col2)), shape=(n_genes, len(barcodes_final)))

    # Remap X1 entries (COO -> filter -> remap) - all bc1 should map
    X1_coo = X1.tocoo()
    r1 = X1_coo.row
    c1 = X1_coo.col
    keep1 = (col_map_X1[c1] != -1)
    if not np.any(keep1):
        logger.warning("No overlapping entries found when aligning X1; returning zero matrix.")
        X1u = sp.csr_matrix((n_genes, len(barcodes_final)))
    else:
        new_row1 = row_map_X1[r1[keep1]]
        new_col1 = col_map_X1[c1[keep1]]
        new_data1 = X1_coo.data[keep1]
        X1u = sp.csr_matrix((new_data1, (new_row1, new_col1)), shape=(n_genes, len(barcodes_final)))

    return X1u.tocsr(), X2u.tocsr(), genes_final, barcodes_final


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
