# PIPseeker Velocity Tools

Tools for generating **RNA velocity–compatible spliced and unspliced matrices** from **PIPseeker** (PIPseq V) outputs.  
This repository includes: A pipeline to build **velocity-ready .h5ad and .loom files**  

---

## Why this exists

Standard RNA velocity methods expect **spliced**, **unspliced** inferred directly from **read-level BAM alignment**.

But for **PIPseq V**, BAM parsing is *incorrect*.

### Why BAM analysis does *not* work for PIPseq V

PIPseeker transforms raw reads into molecule counts using:

- MI correction (binning index + IMI)
- Bin-aware IMI deduplication
- Fragmentation-aware combinatorial inflation modeling
- Multi-mapping resolution

None of this logic survives in the **BAM file**.  
Parsing the BAM yields *read counts*, not *molecule counts*, and treats PIPseq V MIs as UMIs—which is invalid.

### Correct solution: dual-run subtraction

To recover **molecule-level spliced and unspliced** counts:

1. **Run PIPseeker normally** → counts include *exonic + intronic* molecules  
2. **Run PIPseeker with `--exons-only`** on the **raw/unfiltered** matrix → spliced-only molecules  
3. Compute:  
   \
   	ext{unspliced} = 	ext{total} - 	ext{spliced}

This preserves PIPseeker’s counting model and produces valid velocity layers.

This repository automates the process.

---

## Contents

```
pipseeker_velocity_matrices.py     # Build velocity-ready H5AD + LOOM from two PIPseeker runs
```

> `make_pipseeker_test_data.py`  # Generate synthetic PIPseeker test datasets
---

## Build velocity matrices

```bash
python pipseeker_velocity_matrices.py   \
        --total  /path/to/pipseeker_total   \
        --exonic /path/to/pipseeker_exonic_raw   \
        --out-h5ad pipseq_velocity.h5ad   \
        --out-loom pipseq_velocity.loom   \
        --genes-col 0
```

Add scVelo preprocessing:

```bash
--run-scvelo-preproc
```

### Outputs include:

- `spliced` (exonic-only molecules)  
- `unspliced` (total - exonic)  
- `.h5ad` + `.loom` velocity-compatible files

---

## 📘 How it works

1. Load both raw matrices (total & exonic-only)  
2. Align genes and barcodes (union index)  
3. Compute velocity layers  
4. Build an AnnData object with:

```
layers["spliced"]
layers["unspliced"]
```

5. Export to H5AD and LOOM  

---

## ❗ Important Notes

### ✔ Use **RAW/UNFILTERED** exonic matrices
Filtered exons-only matrices may have a different set of barcodes than the total run.

### ✔ PIPseq V MIs are *not* UMIs
They cannot be used for BAM-based UMI counting, exon/intron classification, or velocity reconstruction.



## Contact

For questions, suggestions, or feature requests, contact us at ccrifx@nih.gov
