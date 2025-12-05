# velocity-kit

[![PyPI version](https://badge.fury.io/py/velocity-kit.svg)](https://badge.fury.io/py/velocity-kit)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

**A cross-platform toolkit for building RNA velocity-ready spliced/unspliced matrices** from 10x Genomics, Parse Biosciences, ScaleBio, Fluent BioSciences (PIPseq), and other single-cell technologies.

## Overview

Standard RNA velocity methods expect **spliced** and **unspliced** counts, but many modern single-cell platforms don't directly output these layers. `velocity-kit` provides platform-specific tools to generate velocity-compatible matrices using the **dual-run subtraction method**.

### Supported Platforms

- ✅ **Fluent BioSciences (PIPseq)** - via PIPseeker
- 🚧 **10x Genomics** - Coming soon
- 🚧 **Parse Biosciences** - Coming soon  
- 🚧 **ScaleBio** - Coming soon

## Installation

### From PyPI (recommended)

```bash
pip install velocitykit
```

### From source

```bash
git clone https://github.com/yourusername/velocity-kit.git
cd velocity-kit
pip install -e .
```

### Optional dependencies

To run scVelo preprocessing:

```bash
pip install velocity-kit[scvelo]
```

For development:

```bash
pip install velocity-kit[dev]
```

## Quick Start

### PIPseq (PIPseeker)

```bash
# Generate both h5ad and loom files
velocity-kit prep-pipseq \
  --total /path/to/pipseeker_total_run \
  --exonic /path/to/pipseeker_exons_only_run \
  --out-h5ad output.h5ad \
  --out-loom output.loom

# Or generate only h5ad
velocity-kit prep-pipseq \
  --total /path/to/pipseeker_total_run \
  --exonic /path/to/pipseeker_exons_only_run \
  --out-h5ad output.h5ad

# Or generate only loom
velocity-kit prep-pipseq \
  --total /path/to/pipseeker_total_run \
  --exonic /path/to/pipseeker_exons_only_run \
  --out-loom output.loom
```

## Usage

### Command Structure

```bash
velocity-kit <platform-command> [options]
```

Available platform commands:
- `prep-pipseq` - Prepare velocity matrices from PIPseeker outputs
- `prep-tenx` - Prepare velocity matrices from 10x Genomics outputs (coming soon)
- `prep-parse` - Prepare velocity matrices from Parse Biosciences outputs (coming soon)
- `prep-scalebio` - Prepare velocity matrices from ScaleBio outputs (coming soon)

### PIPseq Detailed Usage

#### Required Arguments

- `--total`: Directory with PIPseeker run that includes introns (total counts)
- `--exonic`: Directory with PIPseeker `--exons-only` run using the RAW/UNFILTERED count matrix
- At least one of:
  - `--out-h5ad`: Output `.h5ad` file path
  - `--out-loom`: Output `.loom` file path

#### Optional Arguments

- `--genes-col`: Column index in `features.tsv` to use as gene ID (default: 0)
- `--run-scvelo-preproc`: Run basic scVelo preprocessing on the AnnData object
- `-v, --verbose`: Increase verbosity level (use `-v` for info, `-vv` for debug)

#### Example

```bash
velocity-kit prep-pipseq \
  --total Analysis/total_run \
  --exonic Analysis/exonic_raw_run \
  --out-h5ad velocity.h5ad \
  --out-loom velocity.loom \
  --run-scvelo-preproc \
  -v
```

### Python API

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

# Align to union of genes and barcodes
X_total_u, X_exon_u, genes_u, bc_u = align_and_union(
    X_total, bc_total, g_total,
    X_exon, bc_exon, g_exon
)

# Build velocity-compatible AnnData
adata = build_velocity_adata(X_total_u, X_exon_u, genes_u, bc_u)

# Save
adata.write_h5ad("output.h5ad")
adata.write_loom("output.loom")
```

## Why Dual-Run Subtraction?

For platforms like PIPseq that use complex molecular counting (MI correction, deduplication, multi-mapping resolution), BAM-based velocity methods are **invalid** because these counting transformations don't survive in the BAM file.

The correct approach:

1. **Run your pipeline normally** → counts include exonic + intronic molecules
2. **Run with exons-only mode** on the **raw/unfiltered** matrix → spliced-only molecules
3. Compute: **unspliced = total - spliced**

This preserves the platform's counting model and produces valid velocity layers.

## Important Notes

⚠️ **For PIPseq**: The `--exonic` directory must point to the **RAW/UNFILTERED** exons-only run.

Do NOT use a filtered exonic matrix, because the called-cell set may not match the total matrix. This will cause barcode mismatches and incorrect velocity estimates.

## Requirements

- Python ≥ 3.7
- anndata ≥ 0.8.0
- h5py ≥ 3.8.0
- numpy ≥ 1.21.0
- pandas ≥ 1.3.0
- scipy ≥ 1.7.0
- tqdm ≥ 4.60.0

Optional:
- scvelo ≥ 0.2.4 (for preprocessing)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

MIT License - see LICENSE file for details.

## Citation

If you use this tool in your research, please cite:

```
[Add citation information here]
```

## Contact

For questions or issues, please open an issue on [GitHub](https://github.com/yourusername/velocity-kit/issues).

## Changelog

### v0.1.0 (Initial Release)
- PIPseq/PIPseeker support
- Modular platform architecture
- Python API for custom workflows
