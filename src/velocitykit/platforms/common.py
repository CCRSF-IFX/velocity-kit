"""Generic dual-run subtraction implementation for velocity-kit.

This module provides a generic implementation that works for any platform
using the dual-run subtraction method (total - exonic = unspliced).
"""

import logging
from pathlib import Path

try:
    from tqdm.auto import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

from ..core import (
    load_10x_mtx,
    align_and_union,
    build_velocity_adata,
    run_scvelo_preprocessing,
)

logger = logging.getLogger(__name__)


def run_dual_subtraction(
    args,
    platform_name: str,
    total_help_text: str,
    exonic_help_text: str,
    subdirectory: str = None
):
    """
    Generic dual-run subtraction pipeline.
    
    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments
    platform_name : str
        Name of the platform for logging (e.g., "PIPseeker", "10x Genomics")
    total_help_text : str
        Description of the total counts input for error messages
    exonic_help_text : str
        Description of the exonic counts input for error messages
    subdirectory : str, optional
        If provided, will look for this subdirectory (e.g., "raw_feature_bc_matrix")
        when the matrix files are not found directly
    """
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
    logger.info(f"NOTE: {total_help_text}")
    logger.info(f"      {exonic_help_text}")

    # Count steps based on outputs
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

    logger.info("Loading total matrix (introns included)...")
    X_total, bc_total, g_total = load_10x_mtx(
        total_dir / "matrix.mtx.gz",
        total_dir / "barcodes.tsv.gz",
        total_dir / "features.tsv.gz",
        genes_col=args.genes_col,
    )
    step_done()

    logger.info("Loading exons-only matrix (RAW/UNFILTERED)...")
    X_exon, bc_exon, g_exon = load_10x_mtx(
        ex_dir / "matrix.mtx.gz",
        ex_dir / "barcodes.tsv.gz",
        ex_dir / "features.tsv.gz",
        genes_col=args.genes_col,
    )
    step_done()

    logger.info("Aligning matrices to union of genes and barcodes...")
    X_total_u, X_exon_u, genes_u, bc_u = align_and_union(
        X_total, bc_total, g_total,
        X_exon, bc_exon, g_exon,
    )
    step_done()

    logger.info("Building velocity-compatible AnnData...")
    adata = build_velocity_adata(X_total_u, X_exon_u, genes_u, bc_u)
    step_done()

    # Write outputs
    if args.out_h5ad:
        out_h5ad = Path(args.out_h5ad)
        logger.info(f"Writing .h5ad to {out_h5ad}")
        adata.write_h5ad(str(out_h5ad))
        step_done()

    if args.out_loom:
        out_loom = Path(args.out_loom)
        logger.info(f"Writing .loom to {out_loom}")
        adata.write_loom(str(out_loom))
        step_done()

    if pbar is not None:
        pbar.close()

    logger.info("✅ Finished building velocity-compatible files.")
    if args.out_h5ad:
        logger.info(f"H5AD: {out_h5ad}")
    if args.out_loom:
        logger.info(f"LOOM: {out_loom}")


def _find_matrix_dir(base_dir: Path, subdirectory: str) -> Path:
    """
    Find the directory containing matrix files.
    
    Looks for matrix.mtx.gz in base_dir first, then in base_dir/subdirectory.
    
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
    matrix_file = base_dir / "matrix.mtx.gz"
    if matrix_file.exists():
        return base_dir
    
    # Try subdirectory
    subdir_matrix = base_dir / subdirectory / "matrix.mtx.gz"
    if subdir_matrix.exists():
        logger.info(f"Found matrix files in {base_dir / subdirectory}")
        return base_dir / subdirectory
    
    raise FileNotFoundError(
        f"Could not find matrix.mtx.gz in {base_dir} or {base_dir / subdirectory}"
    )


def add_standard_arguments(parser, platform_name: str, default_genes_col: int = 0):
    """
    Add standard arguments for dual-run subtraction.
    
    Parameters
    ----------
    parser : argparse.ArgumentParser
        Argument parser to add arguments to
    platform_name : str
        Platform name for help text
    default_genes_col : int
        Default column index for gene IDs in features.tsv
    """
    parser.add_argument(
        "--total",
        required=True,
        help=f"Directory with {platform_name} run that includes introns (total counts).",
    )
    parser.add_argument(
        "--exonic",
        required=True,
        help=(
            f"Directory with {platform_name} exons-only run using the RAW/UNFILTERED "
            "count matrix (before cell calling)."
        ),
    )
    parser.add_argument(
        "--genes-col",
        type=int,
        default=default_genes_col,
        help=f"Column index in features.tsv to use as gene ID (default: {default_genes_col}).",
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
        help="Increase verbosity level (-v, -vv).",
    )
