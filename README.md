# velocity-kit

[![PyPI version](https://badge.fury.io/py/velocity-kit.svg)](https://badge.fury.io/py/velocity-kit)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

Standard RNA velocity methods expect **spliced** and **unspliced** count
layers, but single-cell pipelines expose the information needed to build those
layers in different ways. `velocity-kit` converts platform-native outputs into
velocity-ready H5AD and loom files using the preparation strategy appropriate
for each platform.

`velocity-kit` is not limited to dual-run subtraction. It currently supports both
direct transcript classification and subtraction of matched count runs:

### Supported Platforms

| Platform | Input | Preparation strategy |
| --- | --- | --- |
| **Parse Biosciences** | Split Pipe transcript assignments | Directly classify transcripts with the `exonic` field; no second count run or subtraction |
| **Fluent BioSciences (PIPseq)** | Matched standard and exons-only PIPseeker runs | Dual-run subtraction (`unspliced = total - exonic`) |
| **10x Genomics** | Matched Cell Ranger runs with and without introns | Dual-run subtraction (`unspliced = total - exonic`) |

ScaleBio support is planned but not yet implemented.

## Installation

### From PyPI (recommended)

```bash
pip install velocity-kit
```

### From source

```bash
git clone https://github.com/CCRSF-IFX/velocity-kit.git
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
# Build velocity layers by dual-run subtraction
velocity-kit prep-pipseq \
  --total /path/to/pipseeker_total_run \
  --exonic /path/to/pipseeker_exons_only_run \
  --out-loom output.loom

# Step 2: Generate analysis report (optional)
velocity-kit run-scvelo output.loom -o reports/sample1
```

### 10x Genomics (CellRanger)

```bash
# Build velocity layers by dual-run subtraction
velocity-kit prep-tenx \
  --total /path/to/cellranger_with_introns/raw_feature_bc_matrix \
  --exonic /path/to/cellranger_standard/raw_feature_bc_matrix \
  --out-loom output.loom

# Step 2: Generate analysis report (optional)
velocity-kit run-scvelo output.loom -o reports/sample1
```

### Parse Biosciences (Split Pipe)

Parse data are built directly from transcript assignments; an exons-only
rerun is not required.

```bash
velocity-kit prep-parse \
  --sublibrary /path/to/sublibrary1 \
  --sublibrary /path/to/sublibrary2 \
  --combined-metadata /path/to/combined-output \
  --out-h5ad parse_velocity.h5ad \
  --out-loom parse_velocity.loom
```

`velocity-kit` reads `process/tscp_assignment.csv.gz`, restricts counts to the
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
- `prep-parse` - Prepare velocity matrices directly from Parse Biosciences Split Pipe transcript assignments
- `prep-scalebio` - Prepare velocity matrices from ScaleBio outputs (coming soon)
- `assemble` - Combine multiple loom/H5AD inputs from a sample manifest
- `correct-batch` - Correct batch effects while preserving spliced/unspliced ratios
- `run-scvelo` - Run scVelo analysis and generate a report from loom or H5AD input

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
# Generate velocity-compatible matrices by dual-run subtraction
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
# Method 1: Point to the matched count directories directly
velocity-kit prep-tenx \
  --total cellranger_introns/outs/raw_feature_bc_matrix \
  --exonic cellranger_standard/outs/raw_feature_bc_matrix \
  --out-loom velocity.loom \
  -v

# Generate analysis report
velocity-kit run-scvelo velocity.loom -o reports/sample1

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
provided. `velocity-kit` uses the ordered sublibrary paths recorded in the Split
Pipe log and verifies every local metadata row against the corresponding
combined `__sN` partition. If only the CSV is supplied, a unique exact metadata
match is required. Missing, inconsistent, or ambiguous mappings stop with an
error rather than silently assigning the wrong cell identities. For relocated
results, only each logged path's base folder name is used under the supplied
`--sublibraries-dir`.

### Assemble Multiple Velocity Inputs

Combine multiple `.loom` or `.h5ad` files into one canonical H5AD using a
CSV/TSV manifest. The manifest requires a `path` column; other columns become
per-cell annotations. An optional unique `source` column controls the suffix
used to make cell identifiers unique.

```csv
path,source,sample,batch
E6_rep1.loom,E6_rep1,Cb_E6,run1
E6_rep2.loom,E6_rep2,Cb_E6,run2
E7.loom,E7,Cb_E7,run1
```

```bash
velocity-kit assemble \
  --manifest samples.csv \
  --output combined_velocity.h5ad
```

By default, all inputs must contain exactly the same unique gene identifiers;
gene order may differ and is aligned automatically. Use
`--gene-join intersection` explicitly to retain only genes shared by every
input. `X` is standardized to the `spliced` layer, while both kinetic layers
and assembly provenance are preserved.

### Velocity-Aware Batch Correction

`correct-batch` implements the
[Hansen-ComBat strategy](https://www.hansenlab.org/velocity_batch): it jointly
normalizes total spliced-plus-unspliced abundance, applies ComBat to
`log1p(S + U)`, and reconstructs the two corrected layers using each original
gene/cell spliced fraction. This removes abundance-scale batch effects without
correcting spliced and unspliced independently and distorting their kinetic
relationship.

```bash
velocity-kit correct-batch combined_velocity.h5ad \
  --output combined_velocity_corrected.h5ad \
  --batch-key batch \
  --preserve-key sample cell_type \
  --n-top-genes 2000
```

The batch column must contain at least two batches with at least two cells in
each. Use `--preserve-key` only for biological covariates that should remain in
the model; it cannot include the batch column. Strongly confounded batch and
biological variables cannot be separated reliably. Use `--n-top-genes 0` to
retain all genes, with higher memory use.

The output is a canonical H5AD with normalized, corrected `spliced` and
`unspliced` layers, `X = log1p(spliced)`, pre-correction per-cell QC metrics,
and correction provenance. `run-scvelo` recognizes that provenance and skips a
second normalization and variable-gene-selection pass.

### scVelo Analysis Report

Generate a comprehensive HTML report with QC plots, velocity analysis, and visualizations from a `.loom` or `.h5ad` file.

#### Required Arguments

- `input_path`: Path to an input `.loom` or `.h5ad` file generated by `prep-*` commands

#### Optional Arguments

- `-o, --output-dir`: Output directory for plots and HTML report (default: `scvelo_analysis`)
- `-n, --sample-name`: Sample name for report title (default: derived from input filename)
- `--metadata-file`: Optional cell metadata CSV or TSV
- `--metadata-key`: Unique cell identifier column in the metadata file
- `--adata-key`: Matching column in the loom observations; use `_index` for observation names
- `--color-by COLUMN [COLUMN ...]`: Generate additional UMAPs colored by selected annotation columns
- `--subset-by COLUMN`: Observation/metadata column used to select cells before preprocessing
- `--subset-values VALUE [VALUE ...]`: One or more values retained from `--subset-by`
- `--save-anndata PATH`: Save the fully analyzed AnnData to an `.h5ad` file
- `--analysis-mode {joint,per-group}`: Analyze selected cells together or independently by group
- `--group-by COLUMN`: Column defining independent analyses in `per-group` mode
- `-v, --verbose`: Increase verbosity level (use `-v` for info, `-vv` for debug)

External metadata are joined before the expensive scVelo calculation. Every
input cell must have exactly one matching metadata row. Duplicate identifiers,
missing cells, ambiguous automatic keys, and conflicts with existing loom
annotations stop the run with an error. Extra metadata rows are allowed and
reported. If the identifier column has the same name in both inputs,
`--adata-key` can be omitted.

Cell subsetting is applied after external metadata validation and before
normalization, variable-gene selection, neighborhood construction, or velocity
estimation. `--subset-by` and `--subset-values` must be supplied together, and
every requested value must exist.

`--analysis-mode joint` is the default and creates one neighbor graph and
velocity model across all selected cells. `--analysis-mode per-group` runs the
entire pipeline independently for every `--group-by` value and writes a parent
HTML index linking the group reports. If `--save-anndata analyzed.h5ad` is also
used, per-group files are named `analyzed_<group>.h5ad`; a path containing
`{group}` can be used as an explicit template.

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

# Attach external metadata and add sample/cell-type UMAPs to the report
velocity-kit run-scvelo velocity.loom \
  -o reports/sample1 \
  --metadata-file cell_metadata.csv \
  --metadata-key barcode \
  --color-by sample cell_type

# Analyze only Cb_E6 cells
velocity-kit run-scvelo velocity.loom \
  -o reports/Cb_E6 \
  --metadata-file cell_metadata.csv \
  --metadata-key bc_wells \
  --adata-key bc_wells \
  --subset-by sample \
  --subset-values Cb_E6 \
  --color-by sample

# Analyze the union of multiple samples
velocity-kit run-scvelo velocity.loom \
  -o reports/Cb_E6_E7 \
  --subset-by sample \
  --subset-values Cb_E6 Cb_E7 \
  --color-by sample \
  --save-anndata reports/Cb_E6_E7/analyzed.h5ad

# Run fully independent analyses for every selected sample
velocity-kit run-scvelo velocity.loom \
  -o reports/by_sample \
  --analysis-mode per-group \
  --group-by sample \
  --color-by sample \
  --save-anndata reports/by_sample/analyzed_{group}.h5ad

# Use default output directory and auto-detect sample name
velocity-kit run-scvelo velocity.loom
```

#### Output

The report includes:
- **QC plots**: Input total counts and gene counts calculated before gene filtering,
  plus spliced/unspliced proportions
- **Velocity embeddings**: UMAP with velocity arrows and stream plots
- **Metadata embeddings**: One UMAP for each requested `--color-by` column
- **Clustering**: Leiden community detection (resolution=0.1)
- **Analyzed AnnData**: Optional H5AD containing velocity layers, graphs, embeddings, clusters, and metadata
- **Top velocity genes**: Ranked genes driving velocity patterns
- **HTML report**: All plots combined in an interactive HTML file

#### Features (v0.2.0+)

- **Adaptive parameters**: Automatically adjusts `n_neighbors` and `n_pcs` for small datasets (<50 cells)
- **Robust gene alignment**: Uses scanpy for reliable data loading and processing
- **Leiden clustering**: Modern community detection algorithm (replaces deprecated Louvain)
- **Velocyto compatibility**: Output format matches velocyto standard for downstream tools

### Python API

#### Dual-run subtraction

```python
import scanpy as sc
from velocitykit.platforms.common import build_velocity_adata_from_anndata

# Load matrices using scanpy (recommended in v0.2.0+)
adata_total = sc.read_10x_mtx("total_run/raw_feature_bc_matrix", var_names='gene_symbols', gex_only=True)
adata_exonic = sc.read_10x_mtx("exonic_run/raw_feature_bc_matrix", var_names='gene_symbols', gex_only=True)

# Build velocity-compatible AnnData with proper gene alignment
adata = build_velocity_adata_from_anndata(adata_total, adata_exonic)

# Save in multiple formats
adata.write_h5ad("output.h5ad")
adata.write_loom("output.loom")
```

#### Parse transcript assignments

```python
from pathlib import Path
from velocitykit import build_parse_sublibrary_adata, combine_parse_sublibraries

sublibrary = build_parse_sublibrary_adata(
    Path("sublibrary1/process/tscp_assignment.csv.gz"),
    Path("sublibrary1/all-sample/DGE_filtered/cell_metadata.csv"),
)

adata = combine_parse_sublibraries([sublibrary], labels=["1"])
adata.write_h5ad("parse_velocity.h5ad")
```

For multiple Parse sublibraries, the `prep-parse` CLI is recommended because
it validates the `__sN` mapping against the combined Split Pipe output.

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

# Align the matrices to the total run's genes and filtered barcodes
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

## Preparation Strategies

### Direct Transcript Classification

Parse Split Pipe records one row per assigned transcript and provides an
`exonic` classification. `velocity-kit` uses this information directly:

1. Retain transcripts belonging to the filtered cells in each sublibrary.
2. Assign `exonic=True` transcripts to the `spliced` layer.
3. Assign `exonic=False` transcripts to the `unspliced` layer.
4. Combine sublibraries while reproducing Split Pipe's `__sN` cell identities.

This method uses a single set of Split Pipe results and does not perform
subtraction.

### Dual-Run Subtraction

For platforms that use complex molecular counting (MI correction, deduplication, multi-mapping resolution), BAM-based velocity methods can be **invalid** because these counting transformations don't survive in the BAM file.

The **dual-run subtraction** approach:

1. **Run your pipeline normally** → counts include exonic + intronic molecules
2. **Run with exons-only mode** on the **raw/unfiltered** matrix → spliced-only molecules
3. Compute: **unspliced = total - spliced**

This preserves the platform's counting model and produces valid velocity layers.

### Which Strategy Should I Use?

- **Parse Biosciences**: Use `prep-parse` for direct classification from Split
  Pipe transcript assignments.
- **PIPseq**: Use `prep-pipseq` with matched standard and exons-only runs.
- **10x Genomics**: Use `prep-tenx` with matched Cell Ranger runs with and
  without introns.
- **Other platforms**: A new platform adapter can use direct classifications,
  dual-run subtraction, or another platform-appropriate strategy.

## Important Notes

⚠️ **For PIPseq and 10x dual-run workflows**: The `--exonic` directory
must point to the **RAW/UNFILTERED** exons-only matrix.

Do NOT use a filtered exonic matrix, because the called-cell set may not match the total matrix. This will cause barcode mismatches and incorrect velocity estimates.

⚠️ **For Parse Biosciences**: Prefer passing the combined Split Pipe
output directory to `--combined-metadata`. This allows `velocity-kit` to use the
combine log and validate each sublibrary's `__sN` suffix mapping.

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

If you use this tool in your research, please cite: https://github.com/CCRSF-IFX/velocity-kit

## Contact

For questions or issues, please:
- Email: ccrsfifx@nih.gov
- Open an issue on [GitHub](https://github.com/CCRSF-IFX/velocity-kit/issues)
