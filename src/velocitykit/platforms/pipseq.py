"""PIPseq/PIPseeker platform support for velocity-kit."""

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


def add_arguments(parser):
    """Add PIPseq-specific arguments to the argument parser."""
    parser.add_argument(
        "--total",
        required=True,
        help="Directory with PIPseeker run that includes introns (total counts).",
    )
    parser.add_argument(
        "--exonic",
        required=True,
        help=(
            "Directory with PIPseeker --exons-only run using the RAW/UNFILTERED "
            "count matrix (before cell calling)."
        ),
    )
    parser.add_argument(
        "--genes-col",
        type=int,
        default=0,
        help="Column index in features.tsv to use as gene ID (default: 0).",
    )
    parser.add_argument(
        "--out-h5ad",
        help="Output .h5ad file path (optional if --out-loom is specified).",
    )
    parser.add_argument(
        "--out-loom",
        help="Output .loom file path (optional if --out-h5ad is specified).",
    )
    parser.add_argument(
        "--run-scvelo-preproc",
        action="store_true",
        help="If set, run basic scVelo preprocessing on the AnnData object.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=1,
        help="Increase verbosity level (-v, -vv).",
    )


def run(args):
    """Run the PIPseq velocity matrix preparation pipeline."""
    total_dir = Path(args.total)
    ex_dir = Path(args.exonic)

    # Validate output arguments
    if not args.out_h5ad and not args.out_loom:
        raise ValueError(
            "At least one output format must be specified: --out-h5ad or --out-loom"
        )

    # Directory validation
    if not total_dir.is_dir():
        raise NotADirectoryError(f"--total is not a directory: {total_dir}")
    if not ex_dir.is_dir():
        raise NotADirectoryError(f"--exonic is not a directory: {ex_dir}")

    logger.info("Starting PIPseeker velocity matrix construction.")
    logger.info(
        "NOTE: --exonic should be the RAW/UNFILTERED exons-only count matrix "
        "(before cell calling/filtering)."
    )

    # Count steps based on outputs
    num_outputs = int(bool(args.out_h5ad)) + int(bool(args.out_loom))
    steps = 4 + int(args.run_scvelo_preproc) + num_outputs
    if HAS_TQDM:
        pbar = tqdm(total=steps, desc="Pipeline", ncols=80)
    else:
        pbar = None

    def step_done():
        if pbar is not None:
            pbar.update(1)

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

    if args.run_scvelo_preproc:
        logger.info("Running scVelo preprocessing as requested.")
        run_scvelo_preprocessing(adata)
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
