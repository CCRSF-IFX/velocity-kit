"""PIPseq/PIPseeker platform support for velocity-kit."""

import logging
from .common import add_standard_arguments, run_dual_subtraction

logger = logging.getLogger(__name__)


def add_arguments(parser):
    """Add PIPseq-specific arguments to the argument parser."""
    add_standard_arguments(
        parser,
        platform_name="PIPseeker",
        default_genes_col=0  # PIPseq uses column 0 for gene IDs
    )


def run(args):
    """Run the PIPseq velocity matrix preparation pipeline."""
    run_dual_subtraction(
        args,
        platform_name="PIPseeker",
        total_help_text="--total should be the standard PIPseeker run (includes introns)",
        exonic_help_text="--exonic should be the RAW/UNFILTERED exons-only count matrix (before cell calling/filtering)",
        subdirectory=None  # PIPseq doesn't use subdirectories
    )
