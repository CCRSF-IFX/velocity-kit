"""10x Genomics platform support for velocity-kit."""

import logging
from .common import add_standard_arguments, run_dual_subtraction

logger = logging.getLogger(__name__)


def add_arguments(parser):
    """Add 10x Genomics-specific arguments to the argument parser."""
    add_standard_arguments(
        parser,
        platform_name="10x Genomics CellRanger",
        default_genes_col=1  # 10x uses column 1 for gene symbols
    )


def run(args):
    """Run the 10x Genomics velocity matrix preparation pipeline."""
    run_dual_subtraction(
        args,
        platform_name="10x Genomics",
        total_help_text="--total should be from a CellRanger run with --include-introns flag",
        exonic_help_text="--exonic should be the RAW/UNFILTERED exons-only count matrix (e.g., raw_feature_bc_matrix)",
        subdirectory="raw_feature_bc_matrix"  # 10x uses this subdirectory structure
    )
