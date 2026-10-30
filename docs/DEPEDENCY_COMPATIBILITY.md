# Dependency Compatibility

velocity-kit requires **Python 3.9 or newer**. Python 3.8 is no longer supported.

## Dependency Constraints

The core constraints in `pyproject.toml` and `requirements.txt` are:

| Dependency | Constraint |
|------------|------------|
| anndata | `>=0.8.0,<0.13` |
| h5py | `>=3.10.0` |
| leidenalg | `>=0.8.0` |
| loompy | `>=3.0.6` |
| numpy | `>=1.21.0` |
| pandas | `>=1.3.0` |
| scanpy | `>=1.9.0` |
| scipy | `>=1.7.0` |
| tqdm | `>=4.60.0` |

pip selects dependency versions compatible with the Python interpreter and these
constraints. NumPy 2.x is allowed.

AnnData 0.13 exposes the primary matrix as `layers[None]`, which its Loom
writer passes to Loompy and causes an `AttributeError`. Keep the upper bound
until the Loom writer supports this layer representation.

## Installation

Using Python 3.9 or newer:

```bash
python -m pip install velocity-kit
```

## Testing

GitHub Actions tests Python 3.9, 3.10, 3.11, 3.12, and 3.13 on `ubuntu-24.04`.
Keep the CI matrix, package metadata, and tool targets aligned when changing
Python support. Run CI before releasing dependency changes.
