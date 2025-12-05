"""
pipseq-velocity: Tools for generating RNA velocity-compatible matrices from PIPseeker outputs.

This package provides utilities to convert PIPseeker outputs into spliced/unspliced
count matrices suitable for RNA velocity analysis.
"""

__version__ = "0.1.0"
__author__ = "Shaojun Xie"
__email__ = "xies4@nih.gov"

from .core import (
    load_10x_mtx,
    align_and_union,
    build_velocity_adata,
)

__all__ = [
    "load_10x_mtx",
    "align_and_union",
    "build_velocity_adata",
]
