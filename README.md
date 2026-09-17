# velocity-kit

[![PyPI version](https://badge.fury.io/py/velocity-kit.svg)](https://badge.fury.io/py/velocity-kit)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

Standard RNA velocity methods expect **spliced** and **unspliced** counts, but many modern single-cell platforms don't directly output these layers. `velocity-kit` provides platform-specific tools to generate velocity-compatible matrices using the **dual-run subtraction method**.

### Supported Platforms

- ✅ **Fluent BioSciences (PIPseq)** - via PIPseeker
- ✅ **10x Genomics** - via CellRanger with `--include-introns`
- ✅ **Parse Biosciences** - via Split Pipe transcript assignments

## Installation

### From PyPI (recommended)

```bash
pip install velocitykit
```

### From source

```bash
git clone https://github.com/yourusername/velocitykit.git
cd velocitykit
pip install -e .
```

### Optional dependencies

To run scVelo preprocessing:

```bash
pip install velocitykit[scvelo]
```

For development:

```bash
pip install velocity-kit[dev]
```

## Quick Start

### PIPseq (PIPseeker)

```bash
# Step 1: Generate velocity-compatible matrices
velocity-kit prep-pipseq \
  --total /path/to/pipseeker_total_run \
  --exonic /path/to/pipseeker_exons_only_run \
  --out-loom output.loom

# Step 2: Generate analysis report (optional)
velocity-kit run-scvelo output.loom -o reports/sample1
```

### 10x Genomics (CellRanger)

```bash
# Step 1: Generate velocity-compatible matrices
velocity-kit prep-tenx \
  --total /path/to/cellranger_with_introns/raw_feature_bc_matrix \
  --exonic /path/to/cellranger_standard/raw_feature_bc_matrix \
  --out-loom output.loom

# Step 2: Generate analysis report (optional)
velocity-kit run-scvelo output.loom -o reports/sample1
```

### Parse Biosciences (Split Pipe)

```bash
velocity-kit prep-parse \
  --sublibrary /path/to/sublibrary1 \
  --sublibrary /path/to/sublibrary2 \
  --combined-metadata /path/to/combined-output \
  --out-h5ad parse_velocity.h5ad \
  --out-loom parse_velocity.loom
```

VelocityKit reads `process/tscp_assignment.csv.gz`, restricts counts to the
cells in `all-sample/DGE_filtered/cell_metadata.csv`, and reproduces the
combined output's `__s1`, `__s2`, ... cell identities. When
`--combined-metadata` points to the combined output directory, the Split Pipe
log determines the suffix mapping and the combined metadata validates it. CLI
sublibraries may therefore be supplied in any order.

If sublibrary results have moved since combine mode ran, specify their new
parent location instead. Only the base folder names are taken from the log:

```bash
velocity-kit prep-parse \
  --sublibraries-dir /new/location/of/sublibrary-results \
  --combined-metadata /path/to/combined-output \
  --out-h5ad parse_velocity.h5ad
```

**Note**: You can specify just `--out-h5ad` or just `--out-loom` if you only need one format.

## Usage

### Command Structure

```bash
velocity-kit <platform-command> [options]
```

Available platform commands:
- `prep-pipseq` - Prepare velocity matrices from PIPseeker outputs
- `prep-tenx` - Prepare velocity matrices from 10x Genomics CellRanger outputs
- `prep-parse` - Prepare velocity matrices from Parse Biosciences Split Pipe outputs
- `prep-scalebio` - Prepare velocity matrices from ScaleBio outputs (coming soon)
- `run-scvelo` - Run scVelo analysis and generate comprehensive report from loom file

### PIPseq Detailed Usage

#### Required Arguments

- `--total`: Directory with PIPseeker run that includes introns (total counts)
- `--exonic`: Directory with PIPseeker `--exons-only` run using the RAW/UNFILTERED count matrix
- At least one of:
  - `--out-h5ad`: Output `.h5ad` file path
  - `--out-loom`: Output `.loom` file path

#### Optional Arguments

- `--genes-col`: Column index in `features.tsv` to use as gene ID (default: 0)
- `-v, --verbose`: Increase verbosity level (use `-v` for info, `-vv` for debug)


#### Example

```bash
# Generate velocity-compatible matrices
velocity-kit prep-pipseq \
  --total Analysis/total_run \
  --exonic Analysis/exonic_raw_run \
  --out-loom velocity.loom \
  -v

# Generate comprehensive analysis report
velocity-kit run-scvelo velocity.loom \
  -o reports/sample1 \
  -n Sample1
```

### 10x Genomics Detailed Usage

#### Required Arguments

- `--total`: Directory with CellRanger run using `--include-introns` flag (or path to `raw_feature_bc_matrix`)
- `--exonic`: Directory with standard CellRanger run (exons only). Use RAW/UNFILTERED `raw_feature_bc_matrix`, NOT `filtered_feature_bc_matrix`
- At least one of:
  - `--out-h5ad`: Output `.h5ad` file path
  - `--out-loom`: Output `.loom` file path

#### Optional Arguments

- `--genes-col`: Column index in `features.tsv` to use as gene ID (default: 1 for gene symbols)
- `-v, --verbose`: Increase verbosity level (use `-v` for info, `-vv` for debug)


#### Example

```bash
# Method 1: Point to the count directories directly
velocity-kit prep-tenx \
  --total cellranger_introns/outs/raw_feature_bc_matrix \
  --exonic cellranger_standard/outs/raw_feature_bc_matrix \
  --out-loom velocity.loom \
  -v

# Generate analysis report
velocity-kit analyze velocity.loom -o reports/sample1

# Method 2: Point to the parent directories (will auto-find raw_feature_bc_matrix)
velocity-kit prep-tenx \
  --total cellranger_introns/outs \
  --exonic cellranger_standard/outs \
  --out-loom velocity.loom
```

#### How to Generate the Required CellRanger Runs

1. **Standard run (exonic only)**:
   ```bash
   cellranger count --id=sample_exonic \
     --transcriptome=/path/to/refdata \
     --fastqs=/path/to/fastqs \
     --sample=MySample
   ```

2. **Run with introns**:
   ```bash
   cellranger count --id=sample_with_introns \
     --transcriptome=/path/to/refdata \
     --fastqs=/path/to/fastqs \
     --sample=MySample \
     --include-introns
   ```

### Parse Biosciences Detailed Usage

`prep-parse` follows the Parse Biosciences scVelo workflow: rows marked
`exonic=True` in each transcript-assignment table become the `spliced` layer,
and all other assigned transcripts become the `unspliced` layer. The main
matrix contains their sum. Counts are streamed from the compressed files, so
they do not need to be manually decompressed.

#### Required Arguments

- One input form: repeat `--sublibrary` for every Split Pipe sublibrary output,
  or use `--sublibraries-dir` for a new parent directory containing relocated
  result folders whose base names match those recorded in the combine log.
- `--combined-metadata`: Combined Split Pipe output directory (preferred) or
  its `cell_metadata.csv`.
- At least one of `--out-h5ad` or `--out-loom`.

#### Optional Arguments

- `--gene-column`: Use `gene_name` (default) or stable IDs from `gene`.
- `--chunk-size`: Number of transcript rows processed at once (default:
  1,000,000).

The filtered cell metadata is used instead of manually entered transcript
cutoffs. This reproduces Split Pipe's final called-cell set, including runs
where samples within one sublibrary have different calling thresholds.

Suffixes are never inferred from CLI position when combined results are
provided. VelocityKit uses the ordered sublibrary paths recorded in the Split
Pipe log and verifies every local metadata row against the corresponding
combined `__sN` partition. If only the CSV is supplied, a unique exact metadata
match is required. Missing, inconsistent, or ambiguous mappings stop with an
error rather than silently assigning the wrong cell identities. For relocated
results, only each logged path's base folder name is used under the supplied
`--sublibraries-dir`.

### scVelo Analysis Report

Generate a comprehensive HTML report with QC plots, velocity analysis, and visualizations from a loom file.

#### Required Arguments

- `loom_path`: Path to input `.loom` file (generated by `prep-*` commands)

#### Optional Arguments

- `-o, --output-dir`: Output directory for plots and HTML report (default: `scvelo_analysis`)
- `-n, --sample-name`: Sample name for report title (default: derived from loom filename)
- `-v, --verbose`: Increase verbosity level (use `-v` for info, `-vv` for debug)

#### Requirements

This command requires scvelo and scanpy to be installed:

```bash
pip install scvelo scanpy
# or
pip install velocity-kit[scvelo]
```

#### Example

```bash
# Generate analysis report from loom file
velocity-kit run-scvelo velocity.loom \
  -o reports/sample1 \
  -n Sample1 \
  -v

# Use default output directory and auto-detect sample name
velocity-kit run-scvelo velocity.loom
```

#### Output

The report includes:
- **QC plots**: Total counts, gene counts, spliced/unspliced proportions
- **Velocity embeddings**: UMAP with velocity arrows and stream plots
- **Clustering**: Leiden community detection (resolution=0.1)
- **Top velocity genes**: Ranked genes driving velocity patterns
- **HTML report**: All plots combined in an interactive HTML file

#### Features (v0.2.0+)

- **Adaptive parameters**: Automatically adjusts `n_neighbors` and `n_pcs` for small datasets (<50 cells)
- **Robust gene alignment**: Uses scanpy for reliable data loading and processing
- **Leiden clustering**: Modern community detection algorithm (replaces deprecated Louvain)
- **Velocyto compatibility**: Output format matches velocyto standard for downstream tools

### Python API

```python
from velocitykit import build_velocity_adata_from_anndata
import scanpy as sc

# Load matrices using scanpy (recommended in v0.2.0+)
adata_total = sc.read_10x_mtx("total_run/raw_feature_bc_matrix", var_names='gene_symbols', gex_only=True)
adata_exonic = sc.read_10x_mtx("exonic_run/raw_feature_bc_matrix", var_names='gene_symbols', gex_only=True)

# Build velocity-compatible AnnData with proper gene alignment
adata = build_velocity_adata_from_anndata(adata_total, adata_exonic)

# Save in multiple formats
adata.write_h5ad("output.h5ad")
adata.write_loom("output.loom")
```

#### Legacy API (v0.1.x)

```python
from velocitykit import load_10x_mtx, align_and_union, build_velocity_adata
from pathlib import Path

# Load matrices (manual approach)
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

**Note**: The scanpy-based approach (v0.2.0+) is more robust and handles edge cases better. The legacy API is maintained for backward compatibility.

## Why Dual-Run Subtraction?

For platforms that use complex molecular counting (MI correction, deduplication, multi-mapping resolution), BAM-based velocity methods can be **invalid** because these counting transformations don't survive in the BAM file.

The **dual-run subtraction** approach:

1. **Run your pipeline normally** → counts include exonic + intronic molecules
2. **Run with exons-only mode** on the **raw/unfiltered** matrix → spliced-only molecules
3. Compute: **unspliced = total - spliced**

This preserves the platform's counting model and produces valid velocity layers.

### When to Use Dual-Run Subtraction

- ✅ **PIPseq**: Always use dual-run (BAM-based methods are incorrect)
- ✅ **10x Genomics**: Recommended for consistency, especially with CellRanger ≥7.0
- ⚠️ **Other platforms**: Evaluate whether platform-specific counting differs from simple read counting

## Important Notes

⚠️ **For PIPseq**: The `--exonic` directory must point to the **RAW/UNFILTERED** exons-only run.

Do NOT use a filtered exonic matrix, because the called-cell set may not match the total matrix. This will cause barcode mismatches and incorrect velocity estimates.

## Requirements

### Core Dependencies

- Python ≥ 3.8
- anndata ≥ 0.8.0
- h5py ≥ 3.10.0
- leidenalg ≥ 0.8.0
- loompy ≥ 3.0.6
- numpy ≥ 1.21.0 (supports numpy 2.x)
- pandas ≥ 1.3.0
- scanpy ≥ 1.9.0
- scipy ≥ 1.7.0
- tqdm ≥ 4.60.0

### Optional Dependencies

- scvelo ≥ 0.2.4 (for velocity analysis and reports)

### Compatibility Notes

- ✅ **NumPy 2.x**: Fully supported as of v0.2.0
- ✅ **Pandas 2.x**: Categorical column handling fixed in v0.2.0
- ⚠️ **Python 3.7**: Support dropped in v0.2.0. Use velocity-kit v0.1.x for Python 3.7

## Citation

If you use this tool in your research, please cite: https://github.com/CCRSF-IFX/velocitykit

## Contact

For questions or issues, please:
- Email: ccrsfifx@nih.gov
- Open an issue on [GitHub](https://github.com/CCRSF-IFX/velocity-kit/issues)
