"""Shared readers and validation for velocity-kit inputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import scanpy as sc
import scipy.sparse as sp


SUPPORTED_INPUT_SUFFIXES = {".h5ad", ".loom"}


def _stored_values(matrix) -> np.ndarray:
    """Return explicitly stored matrix values without densifying sparse input."""
    if sp.issparse(matrix):
        return np.asarray(matrix.data)
    return np.asarray(matrix).ravel()


def _validate_abundance_matrix(matrix, name: str) -> None:
    values = _stored_values(matrix)
    if not np.issubdtype(values.dtype, np.number):
        raise ValueError(f"{name} must contain numeric abundance values")
    if not np.isfinite(values).all():
        raise ValueError(f"{name} contains non-finite abundance values")
    if np.any(values < 0):
        raise ValueError(f"{name} contains negative abundance values")


def validate_velocity_adata(adata) -> None:
    """Validate identifiers and abundance matrices shared by all consumers."""
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(
            "Velocity input must contain at least one cell and one gene"
        )
    if not adata.obs_names.is_unique:
        raise ValueError("Velocity input contains duplicate cell identifiers")
    if not adata.var_names.is_unique:
        raise ValueError("Velocity input contains duplicate gene identifiers")

    missing_layers = [
        layer for layer in ("spliced", "unspliced") if layer not in adata.layers
    ]
    if missing_layers:
        raise ValueError(
            "Input is missing required RNA-velocity layers: "
            f"{', '.join(missing_layers)}"
        )
    if adata.X is None:
        raise ValueError("Velocity input is missing its primary X matrix")

    _validate_abundance_matrix(adata.X, "X")
    for layer in ("spliced", "unspliced"):
        _validate_abundance_matrix(adata.layers[layer], f"Layer {layer!r}")


def validate_raw_velocity_counts(adata) -> None:
    """Require an unprocessed count matrix suitable for batch correction."""
    validate_velocity_adata(adata)
    if "velocitykit_batch_correction" in adata.uns:
        raise ValueError(
            "Input already contains velocity-kit batch-correction provenance; "
            "refusing to correct it a second time"
        )

    analyzed_layers = sorted(
        set(adata.layers).intersection({"Ms", "Mu", "velocity"})
    )
    analyzed_embeddings = sorted(
        key for key in adata.obsm if str(key).startswith("X_")
    )
    if analyzed_layers or analyzed_embeddings:
        artifacts = analyzed_layers + analyzed_embeddings
        raise ValueError(
            "Batch correction requires raw counts, but the input appears to be "
            f"analyzed (found: {', '.join(map(str, artifacts))})"
        )

    for name, matrix in (
        ("X", adata.X),
        ("Layer 'spliced'", adata.layers["spliced"]),
        ("Layer 'unspliced'", adata.layers["unspliced"]),
    ):
        values = _stored_values(matrix)
        if not np.allclose(values, np.rint(values), rtol=0, atol=1e-6):
            raise ValueError(
                f"{name} is not an integer-like raw count matrix; batch "
                "correction must start from unnormalized counts"
            )


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
    validate_velocity_adata(adata)
    return adata
