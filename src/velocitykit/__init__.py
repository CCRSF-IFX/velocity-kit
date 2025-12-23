"""
velocity-kit: A cross-platform toolkit for building RNA velocity-ready spliced/unspliced matrices.

Supports multiple single-cell platforms:
- 10x Genomics
- Parse Biosciences
- ScaleBio
- Fluent BioSciences (PIPseq)
- And more...

This package provides utilities to convert platform-specific count matrices into
spliced/unspliced matrices suitable for RNA velocity analysis.
"""

import warnings

# Suppress known deprecation warnings from dependencies
warnings.filterwarnings('ignore', category=UserWarning, module='louvain')
warnings.filterwarnings('ignore', message='pkg_resources is deprecated')

__version__ = "0.1.3"
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
