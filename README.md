# Metagenomics and Viral Genome Analysis Pipeline

## Overview

This repository contains a bioinformatics workflow for processing metagenomic sequencing data, characterising microbial and viral taxonomic composition, recovering selected viral genomes, and conducting phylogenetic analyses.

The workflow is organised into four stages: Phase 1, Phase 1b, Phase 2, and Phase 3. Each stage contains scripts for a defined part of the analysis.

## Pipeline workflow

| Stage | Description | Main outputs |
|---|---|---|
| Phase 1 | Read quality control and taxonomic profiling | Quality-controlled reads, taxonomic classifications and abundance estimates |
| Phase 1b | Taxon-specific read extraction and reference-based genome assembly | Taxon-specific reads, consensus sequences and assembly summaries |
| Phase 2 | Taxonomic composition and abundance analysis | Summary tables, abundance estimates and visualisations |
| Phase 3 | Sequence preparation and phylogenetic analysis | Multiple sequence alignments, phylogenetic trees and tree visualisations |

## Repository structure

```text
metagenomics-analysis-pipeline/
│
├── 01_quality_control_taxonomic_profiling/
│   ├── step_1_QC_bracken_pipeline.sh
│   ├── step_2_annotation_and_combination.py
│   └── README.md
│
├── 01b_taxon_specific_read_extraction_genome_assembly/
│   ├── step_1_Read_extraction_species_with_children.sh
│   ├── step_2_Reference_based_assembly.sh
│   └── README.md
│
├── 02_taxonomic_abundance_analysis/
│   ├── 01_highest_taxon_analysis.R
│   ├── 02_taxonomic_abundance_analysis.R
│   └── README.md
│
├── 03_phylogenetic_analysis/
│   ├── 01_prepare_sequences.sh
│   ├── 02_align_and_trim.sh
│   ├── 03_iqtree_and_reroot.sh
│   ├── 04_extract_alignment_metadata.R
│   ├── 05_ggtree_visualisation_edited.R
│   ├── close_sequences.R
│   └── README.md
│
└── README.md
```

## Phase 1: Read Quality Control and Taxonomic Profiling

This stage processes raw sequencing reads and characterises their taxonomic composition.

### Scripts

**1. `step_1_QC_bracken_pipeline.sh`**

Performs the initial processing of sequencing reads, including:

- Read quality assessment and filtering using Fastp.
- Host read depletion using Bowtie2.
- Taxonomic classification using Kraken2.
- Taxonomic abundance estimation using Bracken.

**2. `step_2_annotation_and_combination.py`**

Processes Bracken species-level outputs by:

- Annotating taxonomic assignments using NCBI taxonomy files.
- Optionally filtering human-associated entries.
- Recalculating relative abundance fractions based on non-host reads.
- Combining individual sample results into a single dataset.

The main combined output is `combined_species_abundance.csv`.

## Phase 1b: Taxon-Specific Read Extraction and Reference-Based Genome Assembly

This stage uses taxonomic classification results to retrieve reads assigned to selected taxa and their descendant taxa, followed by reference-based genome assembly.

### Scripts

**1. `step_1_Read_extraction_species_with_children.sh`**

Extracts taxon-specific reads from paired-end sequencing datasets using KrakenTools. The script:

- Identifies reads assigned to a selected taxon and its descendants.
- Extracts corresponding paired-end reads.
- Checks paired-end input consistency.
- Skips previously completed outputs where applicable.
- Produces an extraction summary (`extraction_summary.csv`).

**2. `step_2_Reference_based_assembly.sh`**

Performs reference-based genome assembly using:

- BWA for read mapping.
- SAMtools for alignment processing.
- iVar for consensus sequence generation.

The workflow maps extracted reads against selected reference sequences and generates consensus FASTA sequences. Consensus generation uses a minimum frequency threshold of 0.5 and a minimum depth of 10.

The script also calculates resolved bases and the percentage of the reference sequence resolved, producing an assembly summary (`assembly_summary_wide.csv`).

## Phase 2: Taxonomic Composition and Abundance Analysis

This stage uses the combined species abundance dataset generated in Phase 1 to examine taxonomic composition and abundance patterns across samples.

### Scripts

**1. `01_highest_taxon_analysis.R`**

Summarises taxonomic abundance at the highest taxonomic level represented in the dataset. It generates:

- Overall abundance summaries.
- Sample-level taxonomic summaries.
- Number of species represented within each highest taxon.
- Doughnut and stacked bar plots.

**2. `02_taxonomic_abundance_analysis.R`**

Provides a more detailed assessment of taxonomic abundance at family or species level. The analysis includes:

- Overall taxonomic abundance summaries.
- Sample-level abundance summaries.
- Identification of the most abundant taxa.
- Configurable selection of taxonomic groups.
- Summary statistics and visualisations.
- A consolidated analysis report.

The taxonomic level and target group can be adjusted through the script's user settings.

## Phase 3: Sequence Preparation and Phylogenetic Analysis

This stage prepares recovered and reference sequences for phylogenetic analysis, constructs maximum-likelihood trees, and generates metadata-annotated visualisations.

### Scripts

**1. `01_prepare_sequences.sh`**

Cleans input nucleotide FASTA files by:

- Removing duplicate sequence headers.
- Standardising nucleotide characters.
- Retaining valid IUPAC ambiguity codes.
- Replacing invalid nucleotide characters with `N`.

**2. `02_align_and_trim.sh`**

Performs sequence alignment and trimming using:

- MAFFT for multiple sequence alignment.
- TrimAl for automated alignment trimming.

The script also performs sequence orientation and alignment checks, including validation of sequence counts and alignment lengths.

**3. `03_iqtree_and_reroot.sh`**

Constructs maximum-likelihood phylogenetic trees using IQ-TREE. The analysis includes:

- Model selection using ModelFinder.
- 1,000 ultrafast bootstrap replicates.
- 1,000 SH-aLRT replicates.
- Rerooting of phylogenetic trees using a specified outgroup.

**4. `04_extract_alignment_metadata.R`**

Matches sequence identifiers in the trimmed alignment to records in a metadata file. It produces matched metadata records and reports sequences with missing or unused metadata.

**5. `05_ggtree_visualisation_edited.R`**

Generates metadata-annotated phylogenetic visualisations using ggtree. Features include:

- Colouring sequence tips according to selected metadata categories.
- Highlighting selected study sequences.
- Annotating relevant clades.
- Customising tree appearance and legends.

**6. `close_sequences.R`**

Identifies sequences that are phylogenetically closest to specified focal sequences using pairwise tree distances. It extracts the closest sequences for further examination.

## Software and dependencies

The following tools and packages are required across the workflow. Requirements depend on the stage being performed.

### Command-line tools

| Tool | Main application |
|---|---|
| Fastp | Read quality control and filtering |
| Bowtie2 | Host read depletion |
| Kraken2 | Taxonomic classification |
| Bracken | Taxonomic abundance estimation |
| KrakenTools | Taxon-specific read extraction |
| BWA | Reference-based read mapping |
| SAMtools | Alignment processing |
| iVar | Consensus sequence generation |
| MAFFT | Multiple sequence alignment |
| TrimAl | Alignment trimming |
| IQ-TREE | Maximum-likelihood phylogenetic inference |
| Newick Utilities | Phylogenetic tree rerooting |

### R packages

Depending on the script, the workflow uses:

- tidyverse
- viridis
- scales
- stringr
- Polychrome
- ape
- ggtree
- ggplot2
- dplyr
- RColorBrewer

Python is required for the taxonomic annotation and combination script.

## Getting started

1. Clone or download the repository.
2. Install the software dependencies required for the stages you intend to run.
3. Review the `README.md` inside each phase directory.
4. Open each script and update the user-defined settings, including input files, output directories, reference sequences, database locations and metadata paths.
5. Run the scripts in the appropriate order.

The general workflow is:

```text
Raw sequencing reads
        |
        v
Phase 1: Quality control and taxonomic profiling
        |
        v
Phase 1b: Taxon-specific read extraction
         and reference-based genome assembly
        |
        v
Phase 2: Taxonomic composition and abundance analysis
        |
        v
Phase 3: Sequence preparation and phylogenetic analysis
```

Phase 1b and Phase 2 use outputs from Phase 1 for different analytical purposes. Phase 1b focuses on recovering selected sequences, whereas Phase 2 examines taxonomic composition and abundance. Phase 3 uses prepared sequences and associated metadata for phylogenetic analysis.

## Reproducibility and configuration

Scripts contain user-defined settings for specifying input files, output directories, reference sequences, taxonomic targets and metadata. These settings should be reviewed before execution.

The workflow is designed to support adaptation to different metagenomic datasets. Reference sequences, taxonomic databases and metadata should correspond to the study being analysed.

## Scope and limitations

Taxonomic classification and read abundance estimates indicate the presence of sequence reads assigned to particular taxa. They do not, on their own, establish active infection or confirm the presence of viable organisms.

Reference-based genome recovery depends on sequencing depth, read quality, reference sequence suitability and genome coverage. Consensus sequences should therefore be interpreted alongside their corresponding coverage and assembly statistics.

## Citation

If you use or adapt this workflow, please cite the repository and acknowledge the software tools used in the relevant analyses. Repository citation details can be added here following publication or release.
