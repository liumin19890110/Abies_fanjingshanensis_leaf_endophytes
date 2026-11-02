# Leaf endophytic fungal communities in *Abies fanjingshanensis*

This repository provides the data tables, representative amplicon sequence variant (ASV) sequences, and original R analysis script associated with a comparison of young (YL) and old (OL) needles of *Abies fanjingshanensis*.

## Repository contents

### Data (`data/`)
- `d1_ASV_table.txt` — ASV read-count table, with 2,760 ASVs as rows and 50 samples as columns.
- `d2_taxonomy_table.txt` — Taxonomic annotations for the 2,760 ASVs.
- `d3_metadata_table.txt` — Sample IDs, needle-age groups (YL/OL), and paired-sample identifiers (25 pairs).
- `d4_traits_table.txt` — Measurements of nine needle chemical and physiological traits for 50 samples.
- `d5_sem_table.txt` — Input data for structural equation modeling, including diversity and community variables.
- `Representative sequences.fasta` — Representative sequences for the 2,760 ASVs.

### Analysis code (`scripts/`)
- `Liu_et_al_Rscript.R` — Original R analysis script, preserved without modifications.

## Running the original script

**The R script is supplied exactly as provided by the authors.** It begins with a machine-specific working-directory setting:

```r
setwd("/Users/liumin/Desktop/R_Input/")
```

It also reads input files by their original bare filenames (for example, `d1_ASV_table.txt`). Consequently, the repository's `data/` and `scripts/` layout is for organization and **the script will not run directly from the repository root without environment setup**. To run it without editing the script, the user must create the working directory specified above (or otherwise arrange the same working directory), place the required input data files directly in that directory, and ensure the required R packages are installed. The original script may also create intermediate outputs that are used by later sections.

For execution on another computer, users may need to adapt the working-directory line and other environment-dependent settings locally; these adaptations are not included in the deposited original script.

## Data integrity

The ASV identifiers in `d1_ASV_table.txt`, `d2_taxonomy_table.txt`, and `Representative sequences.fasta` correspond one-to-one (2,760 identifiers in each file).

## Raw sequencing data

Raw sequencing reads are deposited in NCBI SRA under BioProject accession **PRJNA1522124**.

## Notes

The deposited tables and script are provided for research transparency and reproducibility. See the associated publication for sampling design, analytical methods, and interpretation.
