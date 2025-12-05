# velocity-kit Package Structure

## Overview
Successfully converted the PIPseq-specific velocity tool into a multi-platform package called `velocity-kit`.

## Package Information
- **Package name**: `velocity-kit` (PyPI)
- **Module name**: `velocitykit` (Python import)
- **CLI command**: `velocity-kit`
- **Version**: 0.1.0
- **Author**: Shaojun Xie <xies4@nih.gov>

## Repository Structure

```
scripts/                              # Main repository (Git root)
├── pyproject.toml                    # Modern Python packaging config
├── README.md                         # User documentation
├── requirements.txt                  # Core dependencies
├── src/
│   └── velocitykit/                  # Main package
│       ├── __init__.py               # Package initialization
│       ├── cli.py                    # Multi-platform CLI with subcommands
│       ├── core.py                   # Core velocity matrix functions
│       └── platforms/                # Platform-specific implementations
│           ├── __init__.py
│           ├── common.py             # Shared dual-run subtraction logic
│           ├── pipseq.py             # PIPseq/PIPseeker wrapper (uses common)
│           └── tenx.py               # 10x Genomics wrapper (uses common)
└── test_data/                        # Test datasets
```

## Supported Platforms

### ✅ Implemented
1. **PIPseq (Fluent BioSciences)** - via PIPseeker
   - Subcommand: `velocity-kit prep-pipseq`
   - Default gene column: 0

2. **10x Genomics** - via CellRanger with `--include-introns`
   - Subcommand: `velocity-kit prep-tenx`
   - Default gene column: 1 (gene symbols)
   - Auto-detects `raw_feature_bc_matrix` subdirectories

### 🚧 Planned
3. **Parse Biosciences**
   - Subcommand: `velocity-kit prep-parse` (placeholder)

4. **ScaleBio**
   - Subcommand: `velocity-kit prep-scalebio` (placeholder)

## CLI Architecture

### Main Command Structure
```bash
velocity-kit <platform-command> [options]
```

### Platform Subcommands
Each platform has its own `prep-*` subcommand that handles platform-specific details:
- Argument parsing via `add_arguments(parser)`
- Execution via `run(args)`

### Flexible Output Formats
Users can generate:
- Only `.h5ad`: `--out-h5ad output.h5ad`
- Only `.loom`: `--out-loom output.loom`
- Both: `--out-h5ad output.h5ad --out-loom output.loom`

At least one output format must be specified.

## Key Features

### 1. Consolidated Code Architecture
- **`platforms/common.py`**: Shared implementation for all dual-run subtraction platforms
  - `run_dual_subtraction()`: Generic pipeline that works for any platform
  - `add_standard_arguments()`: Reusable argument parser setup
  - `_find_matrix_dir()`: Smart directory detection (handles subdirectories like `raw_feature_bc_matrix`)
- **Platform wrappers** (`pipseq.py`, `tenx.py`): Thin wrappers that configure platform-specific details
  - Default gene column (0 for PIPseq, 1 for 10x)
  - Help text and logging messages
  - Directory structure (subdirectories for 10x, flat for PIPseq)

### 2. Dual-Run Subtraction Method
- Subtracts exonic counts from total counts
- Preserves platform-specific counting logic
- Avoids invalid BAM-based approaches

### 2. Matrix Alignment
- Automatically aligns matrices to union of genes and barcodes
- Handles mismatches between total and exonic runs
- Warns about suspicious barcode overlaps

### 3. Data Validation
- Checks for negative unspliced counts (clips to 0)
- Validates file existence and dimensions
- Warns about duplicates and empty IDs

### 4. Optional scVelo Preprocessing
- `--run-scvelo-preproc` flag
- Filters genes, normalizes, log-transforms, computes moments
- Prepares data for immediate velocity analysis

## Installation

### For Users
```bash
# From PyPI (when published)
pip install velocity-kit

# From source
git clone <repo-url>
cd velocity-kit/scripts
pip install -e .

# With optional dependencies
pip install velocity-kit[scvelo]  # For preprocessing
pip install velocity-kit[dev]     # For development
```

### For Developers
```bash
cd scripts
pip install -e ".[dev]"
```

## Usage Examples

### PIPseq Example
```bash
velocity-kit prep-pipseq \
  --total Analysis_mm10/Analysis/1_E2 \
  --exonic Analysis_mm10/Analysis/1_E2_exonic_raw \
  --out-h5ad results/1_E2_velocity.h5ad \
  --run-scvelo-preproc \
  -v
```

### 10x Genomics Example
```bash
velocity-kit prep-tenx \
  --total cellranger_introns/outs/raw_feature_bc_matrix \
  --exonic cellranger_standard/outs/raw_feature_bc_matrix \
  --out-h5ad velocity.h5ad \
  --out-loom velocity.loom \
  -v
```

## Python API

```python
from velocitykit import load_10x_mtx, align_and_union, build_velocity_adata
from pathlib import Path

# Load matrices
X_total, bc_total, g_total = load_10x_mtx(
    Path("total_run/matrix.mtx.gz"),
    Path("total_run/barcodes.tsv.gz"),
    Path("total_run/features.tsv.gz")
)

X_exon, bc_exon, g_exon = load_10x_mtx(
    Path("exonic_run/matrix.mtx.gz"),
    Path("exonic_run/barcodes.tsv.gz"),
    Path("exonic_run/features.tsv.gz")
)

# Align and build
X_total_u, X_exon_u, genes_u, bc_u = align_and_union(
    X_total, bc_total, g_total,
    X_exon, bc_exon, g_exon
)

adata = build_velocity_adata(X_total_u, X_exon_u, genes_u, bc_u)

# Save
adata.write_h5ad("output.h5ad")
```

## Dependencies

### Core Requirements (Python ≥3.7)
- anndata ≥ 0.8.0
- h5py ≥ 3.8.0
- numpy ≥ 1.21.0
- pandas ≥ 1.3.0
- scipy ≥ 1.7.0
- tqdm ≥ 4.60.0

### Optional
- scvelo ≥ 0.2.4 (for preprocessing)
- pytest, black, flake8, mypy (for development)

## Files Modified/Created

### New Files
- `scripts/pyproject.toml` - Modern packaging configuration
- `scripts/README.md` - Comprehensive documentation
- `scripts/requirements.txt` - Dependency list
- `scripts/src/velocitykit/__init__.py` - Package init
- `scripts/src/velocitykit/cli.py` - Multi-platform CLI
- `scripts/src/velocitykit/core.py` - Core functions
- `scripts/src/velocitykit/platforms/__init__.py` - Platform module init
- `scripts/src/velocitykit/platforms/pipseq.py` - PIPseq implementation
- `scripts/src/velocitykit/platforms/tenx.py` - 10x Genomics implementation
- `PACKAGE_STRUCTURE.md` - This file

### Removed/Cleaned
- Removed redundant setup.py, setup.cfg from root
- Removed duplicate pyproject.toml, README.md from root
- Cleaned up old pipseq_velocity references

## Next Steps for Publication

### 1. Update Repository URLs
Replace all instances of `yourusername` in:
- pyproject.toml
- README.md

### 2. Test Installation
```bash
cd scripts
pip install -e .
velocity-kit --help
velocity-kit prep-pipseq --help
velocity-kit prep-tenx --help
```

### 3. Build Distribution
```bash
pip install build twine
python -m build
```

### 4. Publish to Test PyPI
```bash
twine upload --repository testpypi dist/*
```

### 5. Publish to PyPI
```bash
twine upload dist/*
```

### 6. Create GitHub Release
- Tag version v0.1.0
- Upload wheels and source distribution
- Write release notes

## Platform Extension Guide

To add a new platform (e.g., Parse Biosciences):

1. **Create module**: `src/velocitykit/platforms/parse.py`
2. **Import common functions**:
   ```python
   from .common import add_standard_arguments, run_dual_subtraction
   ```
3. **Implement `add_arguments()`**:
   ```python
   def add_arguments(parser):
       """Add Parse-specific arguments."""
       add_standard_arguments(
           parser,
           platform_name="Parse Biosciences",
           default_genes_col=0  # or 1, depends on Parse output
       )
   ```
4. **Implement `run()`**:
   ```python
   def run(args):
       """Run the Parse pipeline."""
       run_dual_subtraction(
           args,
           platform_name="Parse Biosciences",
           total_help_text="Description of total input",
           exonic_help_text="Description of exonic input",
           subdirectory=None  # or "some_subdir" if applicable
       )
   ```
5. **Update**: `platforms/__init__.py` to import new module
6. **Update**: `cli.py` to add subcommand and route
7. **Update**: README.md with usage examples
8. **Test**: Create test data and validate

### Benefits of Consolidated Architecture

- ✅ **No code duplication**: ~150 lines of identical code now shared
- ✅ **Easy to add platforms**: New platforms only need ~20 lines of wrapper code
- ✅ **Consistent behavior**: All platforms use exact same pipeline logic
- ✅ **Single point of maintenance**: Bug fixes and improvements benefit all platforms
- ✅ **Type safety**: Shared code means shared testing coverage

## Testing Checklist

- [ ] Install package in clean environment
- [ ] Test `velocity-kit --help`
- [ ] Test `velocity-kit prep-pipseq` with real data
- [ ] Test `velocity-kit prep-tenx` with real data
- [ ] Test with only `--out-h5ad`
- [ ] Test with only `--out-loom`
- [ ] Test with both output formats
- [ ] Test `--run-scvelo-preproc` flag
- [ ] Test verbosity levels (-v, -vv)
- [ ] Verify Python API imports work
- [ ] Check error handling for missing files
- [ ] Validate output file formats

## Known Issues/TODOs

- [ ] Add unit tests
- [ ] Add integration tests
- [ ] Create example datasets
- [ ] Add CI/CD pipeline
- [ ] Generate API documentation
- [ ] Add progress bars for large datasets
- [ ] Consider adding memory profiling
- [ ] Add support for compressed outputs
