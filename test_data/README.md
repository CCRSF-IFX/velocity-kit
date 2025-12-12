

# Test Data Generator for velocity-kit

This script generates synthetic PIPseq-style test data that mimics the dual-run workflow used in RNA velocity analysis.

## Overview

The generated data simulates the realistic scenario where:

- **Total run** (filtered): Contains exonic + intronic reads for **filtered** cell barcodes
- **Exonic run** (raw/unfiltered): Contains exonic-only reads for **all** cell barcodes (including low-quality cells)

This matches the real PIPseq workflow where the exonic-only run uses raw/unfiltered data with more barcodes than the filtered total run.

## Data Generation Logic

### Cell Barcodes

- **Filtered cells** (e.g., 200): High-quality cells that pass filtering
- **Raw cells** (e.g., 500): All cells including low-quality background cells
- The filtered cells are a **subset** of the raw cells (first N barcodes)

### Count Matrices

**Total matrix** (genes × filtered_cells):
- Total counts ~ Poisson(λ_total) per gene-cell
- Includes both exonic and intronic reads

**Exonic matrix** (genes × raw_cells):
- For **filtered cells**: Exonic counts ~ Binomial(total, p_exon)
  - p_exon varies by gene around the specified fraction (e.g., 0.7)
  - Satisfies: exonic ≤ total (elementwise)
- For **background cells**: Exonic counts ~ Poisson(λ_total × p_exon × 0.1)
  - Lower counts representing low-quality cells
  - ~10% of typical exonic counts

## Usage Examples

### Basic Usage (default: 200 filtered, 500 raw cells)

```bash
python make_pipseeker_test_data.py \
    --out-root test_pipseeker_data \
    --n-genes 500 \
    --n-cells-filtered 200 \
    --n-cells-raw 500 \
    --lambda-total 1.2 \
    --lambda-exonic-fraction 0.7 \
    --seed 42
```

### Small Test Dataset

```bash
python make_pipseeker_test_data.py \
    --out-root test_pipseeker_data \
    --n-genes 100 \
    --n-cells-filtered 50 \
    --n-cells-raw 150 \
    --lambda-total 3.0 \
    --seed 42
```

## Testing velocity-kit

Once you've generated test data, run the pipeline:

```bash
velocity-kit prep-pipseq \
    --total test_pipseeker_data/total/ \
    --exonic test_pipseeker_data/exonic_raw/ \
    --out-loom test.loom \
    -vv
```

## Expected Output

- The output will contain **500 cells** (union of all raw barcodes)
- For the **200 filtered cells**: unspliced = total - exonic (positive values)
- For the **300 background cells**: unspliced = 0 (no total counts available)

