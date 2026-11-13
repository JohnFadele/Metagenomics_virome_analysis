# Phase 3: Sequence Preparation, Phylogenetic Analysis and Visualisation

This phase prepares cleaned nucleotide sequences, generates a trimmed multiple sequence alignment, constructs and reroots a maximum-likelihood phylogenetic tree, and produces metadata-linked tree visualisations. It also includes an optional analysis to identify the closest sequences to selected study sequences using phylogenetic distances.

## Scripts and workflow

1. **`01_prepare_sequences.sh`** – Cleans the input FASTA by removing duplicate headers, converting nucleotide characters to uppercase, retaining valid IUPAC ambiguity codes and replacing invalid characters with `N`.

2. **`02_align_and_trim.sh`** – Reorients sequences and performs multiple sequence alignment using MAFFT. It removes MAFFT orientation markers and trims poorly aligned regions using TrimAl. The trimmed FASTA is used for phylogenetic analysis.

3. **`03_iqtree_and_reroot.sh`** – Constructs a maximum-likelihood phylogenetic tree using IQ-TREE, with ModelFinder, 1,000 ultrafast bootstrap replicates and 1,000 SH-aLRT replicates. Both the ML tree and consensus tree are rerooted using the specified outgroup.

4. **`04_extract_alignment_metadata.R`** – Matches sequence identifiers in the trimmed alignment to the `sequence_id` column in the metadata file. It generates matched, missing and unused metadata tables.

5. **`05_ggtree_visualisation_edited.R`** – Produces a metadata-annotated phylogenetic tree. It supports tip colouring by a selected metadata variable, genotype-consistent clade highlighting and visual identification of study samples.

6. **`close_sequences.R`** – Optionally identifies the 20 closest sequences to each specified focal sequence using pairwise phylogenetic distances and generates a combined list of unique closest sequences.

## Requirements

**Command-line tools:** MAFFT, TrimAl, IQ-TREE, Newick Utilities (`nw_reroot`) and standard Unix utilities.

**R packages:** `ape`, `ggtree`, `ggplot2`, `dplyr` and `RColorBrewer`.

## Inputs and execution

The workflow requires an input nucleotide FASTA file and a metadata CSV containing a sequence identifier column (`sequence_id`). Edit the project paths, input filenames, outgroup and other user settings in each script before execution.

Run the scripts in numerical order:

```bash
bash 01_prepare_sequences.sh
bash 02_align_and_trim.sh
bash 03_iqtree_and_reroot.sh

Rscript 04_extract_alignment_metadata.R
Rscript 05_ggtree_visualisation_edited.R
```

Run `close_sequences.R` separately when closest-sequence analysis is required.

## Main outputs

Results are organised into the following project directories:

- `01_sequence_preparation/` – Cleaned FASTA and sequence-cleaning log.
- `02_alignment/` – Oriented, aligned and trimmed FASTA files, TrimAl report and alignment log.
- `03_IQTREE/` – Original and rooted phylogenetic trees and analysis log.
- `04_metadata/` – Alignment-specific metadata, missing and unused metadata reports, and log.
- `05_visualisation/` – Phylogenetic tree figure and closest-sequence result tables.

The rooted maximum-likelihood tree is the primary input for tree visualisation and closest-sequence analysis.
