"""Command-line interface for velocity-kit with platform-specific subcommands."""

import argparse
import sys
import logging
import os

from . import assemble
from .platforms import parsebio, pipseq, tenx

logger = logging.getLogger(__name__)


def get_version():
    """Get package version from metadata."""
    try:
        # Try to import from package
        from . import __version__
        return __version__
    except ImportError:
        # Try modern importlib.metadata (Python 3.8+)
        try:
            from importlib.metadata import version
            return version("velocity-kit")
        except Exception:
            # Fallback to hardcoded version if package not installed
            return "0.3.0"


def setup_logging(verbosity: int = 1):
    """Configure logging based on verbosity level."""
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG

    logging.basicConfig(
        level=level,
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def main():
    """Main entry point for the velocity-kit CLI."""
    parser = argparse.ArgumentParser(
        prog="velocity-kit",
        description=(
            "A cross-platform toolkit for building RNA velocity-ready "
            "spliced/unspliced matrices from 10x Genomics, Parse Biosciences, "
            "ScaleBio, Fluent BioSciences (PIPseq), and other single-cell technologies."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {get_version()}",
    )

    subparsers = parser.add_subparsers(
        title="platforms",
        description="Platform-specific preparation commands",
        dest="platform",
        help="Choose a platform-specific command",
    )

    # PIPseq subcommand
    pipseq_parser = subparsers.add_parser(
        "prep-pipseq",
        help="Prepare velocity matrices from PIPseeker (Fluent BioSciences PIPseq) outputs",
        description=(
            "Build velocity-compatible spliced/unspliced matrices from two "
            "PIPseeker runs.\n\n"
            "IMPORTANT:\n"
            "  --total  must point to the run that includes introns "
            "           (total = exonic + intronic).\n"
            "  --exonic must point to the RAW/UNFILTERED exons-only run.\n"
            "           Do NOT use a filtered exonic matrix, because the\n"
            "           called-cell set may not match the total matrix."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    pipseq.add_arguments(pipseq_parser)

    # 10x Genomics subcommand
    tenx_parser = subparsers.add_parser(
        "prep-tenx",
        help="Prepare velocity matrices from 10x Genomics CellRanger outputs",
        description=(
            "Build velocity-compatible spliced/unspliced matrices from two "
            "10x Genomics CellRanger runs.\n\n"
            "IMPORTANT:\n"
            "  --total  must point to a CellRanger run with --include-introns flag\n"
            "           (total = exonic + intronic).\n"
            "  --exonic must point to the RAW/UNFILTERED exons-only count matrix\n"
            "           (e.g., raw_feature_bc_matrix directory).\n"
            "           Do NOT use filtered_feature_bc_matrix."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    tenx.add_arguments(tenx_parser)

    # Parse Biosciences subcommand
    parse_parser = subparsers.add_parser(
        "prep-parse",
        help="Prepare velocity matrices from Parse Biosciences Split Pipe outputs",
        description=(
            "Build velocity-compatible spliced/unspliced matrices directly from "
            "Split Pipe transcript assignments. The combine log and combined "
            "cell metadata determine and validate each sublibrary's __sN suffix."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parsebio.add_arguments(parse_parser)

    # ScaleBio subcommand (placeholder)
    scalebio_parser = subparsers.add_parser(
        "prep-scalebio",
        help="Prepare velocity matrices from ScaleBio outputs (Coming soon)",
    )
    scalebio_parser.add_argument("--placeholder", help="Coming soon")

    # Multi-input assembly subcommand
    assemble_parser = subparsers.add_parser(
        "assemble",
        help="Assemble multiple loom/H5AD inputs into one canonical H5AD",
        description=(
            "Read input paths and per-input annotations from a CSV/TSV manifest, "
            "validate spliced/unspliced layers and gene identities, make cell IDs "
            "unique by source, and write one combined H5AD."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    assemble.add_arguments(assemble_parser)

    # scVelo Analysis subcommand
    scvelo_parser = subparsers.add_parser(
        "run-scvelo",
        help="Run scVelo analysis and generate a report from loom or H5AD input",
        description=(
            "Run a complete scVelo pipeline on a .loom or .h5ad file and generate an HTML "
            "report with QC plots, velocity analysis, and visualizations.\n\n"
            "This command requires scvelo to be installed:\n"
            "  pip install scvelo\n"
            "  # or\n"
            "  pip install velocity-kit[scvelo]"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    scvelo_parser.add_argument(
        "input_path",
        help="Path to an input .loom or .h5ad file generated by prep-* commands"
    )
    scvelo_parser.add_argument(
        "-o", "--output-dir",
        default="scvelo_analysis",
        help="Output directory for plots and HTML report (default: scvelo_analysis)"
    )
    scvelo_parser.add_argument(
        "-n", "--sample-name",
        default=None,
        help="Sample name for report title (default: derived from input filename)"
    )
    scvelo_parser.add_argument(
        "--metadata-file",
        default=None,
        help=(
            "Optional cell metadata CSV or TSV. Every input cell must match one "
            "unique metadata row."
        ),
    )
    scvelo_parser.add_argument(
        "--metadata-key",
        default=None,
        help=(
            "Unique cell identifier column in --metadata-file. If omitted, "
            "velocity-kit requires one unambiguous shared identifier column."
        ),
    )
    scvelo_parser.add_argument(
        "--adata-key",
        default=None,
        help=(
            "Matching cell identifier column in the input AnnData observations. "
            "Use '_index' to match AnnData observation names."
        ),
    )
    scvelo_parser.add_argument(
        "--color-by",
        nargs="+",
        default=None,
        metavar="COLUMN",
        help=(
            "One or more observation/metadata columns used to create additional "
            "UMAP plots (for example: --color-by sample cell_type)."
        ),
    )
    scvelo_parser.add_argument(
        "--subset-by",
        default=None,
        metavar="COLUMN",
        help=(
            "Observation/metadata column used to subset cells before scVelo "
            "preprocessing (for example: --subset-by sample)."
        ),
    )
    scvelo_parser.add_argument(
        "--subset-values",
        nargs="+",
        default=None,
        metavar="VALUE",
        help=(
            "One or more values retained from --subset-by "
            "(for example: --subset-values Cb_E6)."
        ),
    )
    scvelo_parser.add_argument(
        "--save-anndata",
        default=None,
        metavar="PATH",
        help=(
            "Optional .h5ad path for the fully analyzed AnnData, including "
            "velocity layers, graphs, embeddings, clusters, and metadata."
        ),
    )
    scvelo_parser.add_argument(
        "--analysis-mode",
        choices=("joint", "per-group"),
        default="joint",
        help=(
            "Analyze all selected cells together (joint, default) or run a "
            "separate complete analysis for every --group-by value."
        ),
    )
    scvelo_parser.add_argument(
        "--group-by",
        default=None,
        metavar="COLUMN",
        help="Observation/metadata column used with --analysis-mode per-group",
    )
    scvelo_parser.add_argument(
        "-v", "--verbose",
        action="count",
        default=1,
        help="Increase verbosity (-v for INFO, -vv for DEBUG)"
    )

    # Parse arguments
    args = parser.parse_args()

    if not args.platform:
        parser.print_help()
        sys.exit(1)

    # Setup logging
    verbosity = getattr(args, "verbose", 1)
    setup_logging(verbosity)

    # Route to platform-specific handler
    if args.platform == "prep-pipseq":
        pipseq.run(args)
    elif args.platform == "prep-tenx":
        tenx.run(args)
    elif args.platform == "prep-parse":
        parsebio.run(args)
    elif args.platform == "prep-scalebio":
        logger.error("ScaleBio support coming soon!")
        sys.exit(1)
    elif args.platform == "assemble":
        assemble.run(args)
    elif args.platform == "run-scvelo":
        run_scvelo(args)
    else:
        parser.print_help()
        sys.exit(1)


def run_scvelo(args):
    """Run scVelo analysis and report generation with optional dependency handling."""
    try:
        from .scvelo_report import run_scvelo_and_generate_report
    except ImportError as e:
        logger.error(
            "scVelo dependencies not found. Please install with:\n"
            "  pip install scvelo scanpy\n"
            "  # or\n"
            "  pip install velocity-kit[scvelo]"
        )
        logger.debug(f"Import error: {e}")
        sys.exit(1)

    # Validate input file exists
    if not os.path.exists(args.input_path):
        logger.error(f"Input file not found: {args.input_path}")
        sys.exit(1)

    if (args.metadata_key or args.adata_key) and not args.metadata_file:
        logger.error("--metadata-key and --adata-key require --metadata-file")
        sys.exit(1)

    if args.metadata_file and not os.path.isfile(args.metadata_file):
        logger.error(f"Metadata file not found: {args.metadata_file}")
        sys.exit(1)

    if bool(args.subset_by) != bool(args.subset_values):
        logger.error("--subset-by and --subset-values must be provided together")
        sys.exit(1)

    if args.save_anndata and not args.save_anndata.lower().endswith(".h5ad"):
        logger.error("--save-anndata must use the .h5ad extension")
        sys.exit(1)

    if args.analysis_mode == "per-group" and not args.group_by:
        logger.error("--group-by is required with --analysis-mode per-group")
        sys.exit(1)
    if args.analysis_mode == "joint" and args.group_by:
        logger.error("--group-by is only valid with --analysis-mode per-group")
        sys.exit(1)

    input_suffix = os.path.splitext(args.input_path)[1].lower()
    if input_suffix not in {".loom", ".h5ad"}:
        logger.error(
            f"Unsupported input format {input_suffix or '<none>'!r}; "
            "expected .loom or .h5ad"
        )
        sys.exit(1)

    logger.info(f"Generating scVelo analysis report for: {args.input_path}")
    logger.info(f"Output directory: {args.output_dir}")

    try:
        report_path = run_scvelo_and_generate_report(
            input_path=args.input_path,
            output_dir=args.output_dir,
            sample_name=args.sample_name,
            metadata_file=args.metadata_file,
            metadata_key=args.metadata_key,
            adata_key=args.adata_key,
            color_by=args.color_by,
            subset_by=args.subset_by,
            subset_values=args.subset_values,
            save_anndata=args.save_anndata,
            analysis_mode=args.analysis_mode,
            group_by=args.group_by,
        )
        logger.info(f"✓ Analysis report successfully generated: {report_path}")
    except Exception as e:
        logger.error(f"Failed to generate analysis report: {e}")
        logger.debug("Full traceback:", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
