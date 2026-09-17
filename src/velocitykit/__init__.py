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

__version__ = "0.3.0"
__author__ = "Shaojun Xie"
__email__ = "xies4@nih.gov"

from .core import (
    load_10x_mtx,
    align_and_union,
    build_velocity_adata,
)
from .platforms.parsebio import (
    build_parse_sublibrary_adata,
    combine_parse_sublibraries,
)

__all__ = [
    "load_10x_mtx",
    "align_and_union",
    "build_velocity_adata",
    "build_parse_sublibrary_adata",
    "combine_parse_sublibraries",
]
