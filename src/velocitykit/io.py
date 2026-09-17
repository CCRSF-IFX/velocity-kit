"""Shared readers and validation for velocity-kit inputs."""

from __future__ import annotations

from pathlib import Path

import scanpy as sc


SUPPORTED_INPUT_SUFFIXES = {".h5ad", ".loom"}


def read_velocity_input(input_path: str):
    """Read a supported velocity AnnData input and validate kinetic layers."""
    path = Path(input_path)
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    if path.suffix.lower() not in SUPPORTED_INPUT_SUFFIXES:
        supported = ", ".join(sorted(SUPPORTED_INPUT_SUFFIXES))
        raise ValueError(
            f"Unsupported input format {path.suffix or '<none>'!r}; expected {supported}"
        )

    adata = sc.read(path)
    missing_layers = [
        layer for layer in ("spliced", "unspliced") if layer not in adata.layers
    ]
    if missing_layers:
        raise ValueError(
            "Input is missing required RNA-velocity layers: "
            f"{', '.join(missing_layers)}"
        )
    return adata
