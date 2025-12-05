"""Command-line interface for velocity-kit with platform-specific subcommands."""

import argparse
import sys
import logging

from .platforms import pipseq

logger = logging.getLogger(__name__)


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
        version="%(prog)s 0.1.0",
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

    # 10x Genomics subcommand (placeholder)
    tenx_parser = subparsers.add_parser(
        "prep-tenx",
        help="Prepare velocity matrices from 10x Genomics outputs (Coming soon)",
    )
    tenx_parser.add_argument("--placeholder", help="Coming soon")

    # Parse Biosciences subcommand (placeholder)
    parse_parser = subparsers.add_parser(
        "prep-parse",
        help="Prepare velocity matrices from Parse Biosciences outputs (Coming soon)",
    )
    parse_parser.add_argument("--placeholder", help="Coming soon")

    # ScaleBio subcommand (placeholder)
    scalebio_parser = subparsers.add_parser(
        "prep-scalebio",
        help="Prepare velocity matrices from ScaleBio outputs (Coming soon)",
    )
    scalebio_parser.add_argument("--placeholder", help="Coming soon")

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
        logger.error("10x Genomics support coming soon!")
        sys.exit(1)
    elif args.platform == "prep-parse":
        logger.error("Parse Biosciences support coming soon!")
        sys.exit(1)
    elif args.platform == "prep-scalebio":
        logger.error("ScaleBio support coming soon!")
        sys.exit(1)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
