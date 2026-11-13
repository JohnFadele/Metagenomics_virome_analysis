# Metagenomics and Viral Genome Analysis Pipeline

## Overview

This repository contains a multi-stage bioinformatics workflow for processing metagenomic sequencing data, characterising microbial taxonomic composition, recovering target viral sequences and conducting phylogenetic analysis.

The pipeline combines quality control, host-read depletion, taxonomic classification, abundance analysis, taxon-specific read extraction, reference-based genome assembly and phylogenetic visualisation. It is organised into three phases, with individual scripts and configuration settings provided for each stage.

## Pipeline workflow

### Phase 1: Read Quality Control and Taxonomic Profiling

This phase processes raw paired-end sequencing reads and generates taxonomic abundance estimates.

The workflow includes:

- Read quality control and filtering using **fastp**.
- Host-read depletion using **Bowtie2**.
- Taxonomic classification using **Kraken2**.
- Abundance estimation using **Bracken**.
- Taxonomic annotation using NCBI taxonomy files (`nodes.dmp` and `names.dmp`).
- Optional removal of human-associated entries and recalculation of relative abundance.
- Consolidation of annotated sample-level results into a combined abundance dataset.

**Main output:** `combined_species_abundance.csv`

### Phase 2: Taxonomic Composition and Abundance Analysis

This phase analyses the consolidated abundance dataset to describe taxonomic composition across samples.

The R scripts perform:

- Overall and sample-level summaries of major taxonomic groups.
- Visualisation of taxonomic composition using doughnut and stacked-bar charts.
- Configurable abundance analysis at family or species level.
- Comparison of overall and sample-level taxonomic abundance.
- Generation of summary tables, publication-quality figures and an analysis report.

**Main outputs:** Taxonomic summary tables, figures and `analysis_report.csv`.

### Phase 3: Sequence Preparation and Phylogenetic Analysis

This phase prepares selected nucleotide sequences for phylogenetic analysis and generates metadata-linked visualisations.

The workflow includes:

- FASTA cleaning, duplicate-header removal and nucleotide character validation.
- Sequence reorientation and multiple sequence alignment using **MAFFT**.
- Alignment trimming using **TrimAl**.
- Maximum-likelihood phylogenetic tree construction using **IQ-TREE**.
- Tree rerooting using **Newick Utilities**.
- Extraction and matching of alignment-specific metadata.
- Phylogenetic tree visualisation using **ggtree**.
- Optional identification of the closest sequences to selected focal sequences using phylogenetic distances.

**Main outputs:** Cleaned and trimmed FASTA files, phylogenetic trees, metadata tables, tree visualisations and closest-sequence tables.

## Repository structure

```text
metagenomics-analysis-pipeline/
│
├── 01_quality_control_taxonomic_profiling/
│   ├── step_1_QC_bracken_pipeline.sh
│   ├── step_2_annotation_and_combination.py
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

## Requirements

The pipeline uses a combination of command-line bioinformatics tools, Python and R.

**Command-line tools**
- fastp
- Bowtie2
- Kraken2
- Bracken
- KrakenTools
- BWA
- SAMtools
- iVar
- MAFFT
- TrimAl
- IQ-TREE
- Newick Utilities

**Programming environments**
- Python 3
- R

**R packages**
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

Additional Python dependencies and tool-specific requirements should be installed according to the scripts and software documentation.

## Getting started

1. Clone or download this repository.
2. Install the required software and dependencies.
3. Prepare the input FASTQ files, reference databases, taxonomy files, reference FASTA sequences and metadata as required by each phase.
4. Edit the user configuration sections in the scripts to specify the correct input files, project directories, database locations and output directories.
5. Run the scripts in the appropriate order, beginning with Phase 1.

Each phase has its own README file containing more detailed instructions, input requirements and expected outputs.

## Reproducibility and configuration

The scripts include user-defined configuration sections for specifying project directories, input files, output locations, reference sequences, taxonomic targets and analysis parameters. These settings should be reviewed and updated before execution.

Intermediate files, summary tables, figures and log files are retained where specified by the individual scripts to support quality assessment, troubleshooting and reproducibility.

## Scope and limitations

The workflow supports metagenomic taxonomic profiling and downstream sequence-based genomic analysis. Taxonomic read assignments and abundance estimates should be interpreted as computational classifications; they do not, on their own, establish active infection or disease.

Genome recovery and phylogenetic analysis depend on sequencing depth, read quality, reference availability and sequence coverage. Results should therefore be assessed using appropriate quality-control and validation criteria.

## Citation

If you use or adapt this pipeline in your research, please cite this repository and acknowledge the bioinformatics tools and databases used in your analysis.

**Repository:** [Insert GitHub repository URL]# Metagenomics_virome_analysis
