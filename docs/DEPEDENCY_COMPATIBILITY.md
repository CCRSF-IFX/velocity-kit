# Dependency Compatibility Matrix

This document provides detailed compatibility information for velocity-kit dependencies across Python versions.

## Summary

velocity-kit supports **Python 3.7-3.12** with carefully selected dependency versions.

## Detailed Compatibility Matrix

### h5py

h5py is used for reading/writing HDF5 files (`.h5ad` format).

| h5py Version | Python 3.7 | Python 3.8 | Python 3.9 | Python 3.10 | Python 3.11 | Python 3.12 | Notes |
|--------------|------------|------------|------------|-------------|-------------|-------------|-------|
| 3.0.x - 3.7.x | ✅         | ✅         | ✅         | ✅          | ✅          | ❌          | Stable, widely used |
| 3.8.x - 3.10.x | ❌        | ✅         | ✅         | ✅          | ✅          | ✅          | Requires Python ≥3.8 |
| 3.11.x+       | ❌        | ✅         | ✅         | ✅          | ✅          | ✅          | Latest versions |

**Constraint with environment markers**:
```toml
"h5py>=3.0.0,<3.8.0; python_version < '3.8'",  # Python 3.7 gets older h5py
"h5py>=3.10.0; python_version >= '3.8'",       # Python 3.8+ gets latest
```

- ✅ Python 3.7 gets h5py 3.0-3.7 (maximum compatibility)
- ✅ Python 3.8-3.12 get h5py 3.10+ (latest features, best performance)
- ✅ No compromise - each Python version gets the best h5py it can use

**Previous approach** (single constraint): `h5py>=3.0.0,<3.8.0`
- ✅ Works on Python 3.7-3.12
- ❌ Python 3.12 users stuck with older h5py 3.7.x
- ⚠️ Users may see "newer version available" warnings

### NumPy

NumPy is the core numerical library used by all scientific Python packages.

| NumPy Version | Python 3.7 | Python 3.8 | Python 3.9 | Python 3.10 | Python 3.11 | Python 3.12 | Notes |
|---------------|------------|------------|------------|-------------|-------------|-------------|-------|
| 1.21.x        | ✅         | ✅         | ✅         | ✅          | ✅          | ❌          | Last version with universal wheels |
| 1.22.x - 1.26.x | ✅       | ✅         | ✅         | ✅          | ✅          | ✅          | Stable, pre-2.0 API |
| 2.0.x - 2.2.x | ❌         | ❌         | ❌         | ✅          | ✅          | ✅          | **Major API changes** |

**Constraint**: `numpy>=1.21.0,<2.0.0`
- ✅ Works on Python 3.7-3.12
- ✅ Avoids NumPy 2.x breaking changes
- ✅ Allows 1.21.x - 1.26.x (5+ years of stable releases)

**Why exclude NumPy 2.x?**
- NumPy 2.0 dropped Python 3.7-3.9 support
- Major API changes that may break existing code
- Not widely adopted yet (as of Dec 2024)

### loompy

loompy is used for writing `.loom` format files for velocity analysis.

| loompy Version | Python 3.7 | Python 3.8 | Python 3.9 | Python 3.10 | Python 3.11 | Python 3.12 | Notes |
|----------------|------------|------------|------------|-------------|-------------|-------------|-------|
| 3.0.6+         | ✅         | ✅         | ✅         | ✅          | ✅          | ✅          | Stable |

**Constraint**: `loompy>=3.0.6`
- ✅ Works on all Python versions 3.7-3.12
- No upper bound needed (stable API)

### pandas

pandas is used for handling barcode and gene metadata.

| pandas Version | Python 3.7 | Python 3.8 | Python 3.9 | Python 3.10 | Python 3.11 | Python 3.12 | Notes |
|----------------|------------|------------|------------|-------------|-------------|-------------|-------|
| 1.3.x          | ✅         | ✅         | ✅         | ✅          | ❌          | ❌          | Minimum for Python 3.7 |
| 1.4.x - 1.5.x  | ✅         | ✅         | ✅         | ✅          | ✅          | ❌          | |
| 2.0.x+         | ❌         | ✅         | ✅         | ✅          | ✅          | ✅          | Dropped Python 3.7 |

**Constraint**: `pandas>=1.3.0`
- ✅ Python 3.7 gets pandas 1.3.x - 1.5.x
- ✅ Python 3.8+ can get pandas 2.x
- No upper bound needed (pip resolver handles it)

### scipy

scipy is used for sparse matrix operations.

| scipy Version | Python 3.7 | Python 3.8 | Python 3.9 | Python 3.10 | Python 3.11 | Python 3.12 | Notes |
|---------------|------------|------------|------------|-------------|-------------|-------------|-------|
| 1.7.x         | ✅         | ✅         | ✅         | ✅          | ✅          | ❌          | Minimum for Python 3.7 |
| 1.8.x - 1.11.x | ✅        | ✅         | ✅         | ✅          | ✅          | ✅          | |
| 1.12.x+       | ❌         | ✅         | ✅         | ✅          | ✅          | ✅          | Dropped Python 3.7 |

**Constraint**: `scipy>=1.7.0`
- ✅ Python 3.7 gets scipy 1.7.x - 1.11.x
- ✅ Python 3.8+ can get latest versions
- No upper bound needed (pip resolver handles it)

### anndata

anndata is the core data structure for single-cell analysis.

| anndata Version | Python 3.7 | Python 3.8 | Python 3.9 | Python 3.10 | Python 3.11 | Python 3.12 | Notes |
|-----------------|------------|------------|------------|-------------|-------------|-------------|-------|
| 0.8.x           | ✅         | ✅         | ✅         | ✅          | ✅          | ✅          | Stable API |
| 0.9.x - 0.10.x  | ✅         | ✅         | ✅         | ✅          | ✅          | ✅          | |

**Constraint**: `anndata>=0.8.0`
- ✅ Works on all Python versions 3.7-3.12
- No upper bound needed (stable API)

### tqdm

tqdm provides progress bars for long-running operations.

| tqdm Version | Python 3.7 | Python 3.8 | Python 3.9 | Python 3.10 | Python 3.11 | Python 3.12 | Notes |
|--------------|------------|------------|------------|-------------|-------------|-------------|-------|
| 4.60.0+      | ✅         | ✅         | ✅         | ✅          | ✅          | ✅          | Universal support |

**Constraint**: `tqdm>=4.60.0`
- ✅ Works on all Python versions 3.7-3.12
- No upper bound needed (stable API)

## Recommended Installation

### For Python 3.7 users:
```bash
pip install velocity-kit
# This will install:
# - h5py~=3.7.0
# - numpy~=1.21.6
# - pandas~=1.5.3
# - scipy~=1.7.3
```

### For Python 3.8-3.12 users:
```bash
pip install velocity-kit
# This will install:
# - h5py>=3.8.0 (or 3.7.x with our constraint)
# - numpy~=1.26.x (pre-2.0)
# - pandas>=2.0.0
# - scipy>=1.12.0
```

## Version Strategy

Our version constraints follow this strategy:

1. **Lower bound**: Minimum version that provides required features
2. **Upper bound**: Only where necessary to avoid breaking changes
3. **Let pip resolve**: For most packages, let pip's dependency resolver pick the best version

### When to add upper bounds:

- ✅ **h5py<3.8.0**: Required for Python 3.7 support
- ✅ **numpy<2.0.0**: Avoid major API changes and dropped Python support
- ❌ **pandas**: Let pip resolver handle it (pandas 2.x is compatible)
- ❌ **scipy**: Let pip resolver handle it (API is stable)
- ❌ **anndata**: Let pip resolver handle it (API is stable)

## Testing Strategy

Our GitHub Actions CI tests on:
- Python 3.7 (ubuntu-20.04)
- Python 3.8-3.12 (ubuntu-latest)

This ensures compatibility across all supported Python versions.

## Updating Dependencies

When considering dependency updates:

1. **Check Python support**: Use [PyPI](https://pypi.org) to verify which Python versions are supported
2. **Test locally**: Use `tox` or `nox` to test on multiple Python versions
3. **Update constraints**: Adjust `pyproject.toml` if needed
4. **Update this document**: Keep the compatibility matrix current
5. **Run CI**: GitHub Actions will test all Python versions

## Future Considerations

### When to drop Python 3.7 support?

Python 3.7 reached end-of-life on **June 27, 2023**. Consider dropping support when:

- ✅ Major dependencies (NumPy 2.x, pandas 3.x) require newer Python
- ✅ User base has migrated to Python 3.8+
- ✅ New features require Python 3.8+ syntax (e.g., walrus operator `:=`)

**Benefits of dropping Python 3.7**:
- Use h5py ≥3.8.0 (better performance, bug fixes)
- Use NumPy 2.x (if needed)
- Simplified dependency constraints
- Use Python 3.8+ features (walrus operator, positional-only parameters)

**Current status**: We maintain Python 3.7 support because:
- Users may be on older HPC systems
- All features work fine on Python 3.7
- No significant maintenance burden
