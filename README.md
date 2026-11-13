# Metagenomics and Viral Genome Analysis Pipeline

## Overview

This repository contains a bioinformatics workflow for processing metagenomic sequencing data, characterising taxonomic composition, recovering selected microbial and viral genomes, and performing phylogenetic analysis.

The workflow is organised into four steps, covering read quality control, taxonomic profiling, abundance analysis, taxon-specific read extraction, reference-based genome recovery, and phylogenetic analysis.

## Pipeline Workflow

The pipeline consists of four main steps:

| Step | Description |
|---|---|
| Phase 1 | Read Quality Control and Taxonomic Profiling |
| Phase 1b | Taxonomic Composition and Abundance Analysis |
| Phase 2 | Taxon-Specific Read Extraction and Reference-Based Genome Assembly |
| Phase 3 | Sequence Preparation and Phylogenetic Analysis |

### Phase 1: Read Quality Control and Taxonomic Profiling

This phase processes raw sequencing reads and determines the taxonomic composition of each sample.

**Main operations:**
- Quality filtering of raw reads using fastp.
- Host-read depletion using Bowtie2.
- Taxonomic classification using Kraken2.
- Taxonomic abundance estimation using Bracken.
- Annotation of Bracken results using NCBI taxonomy files.
- Optional filtering of *Homo sapiens* entries.
- Recalculation of relative abundance based on non-host reads.
- Combination of annotated species-level abundance results.

**Scripts:**
- `step_1_QC_bracken_pipeline.sh`
- `step_2_annotation_and_combination.py`

**Main output:**
- `combined_species_abundance.csv`

### Phase 1b: Taxonomic Composition and Abundance Analysis

This phase analyses the combined abundance dataset generated in Phase 1 to describe and visualise the taxonomic composition of the samples.

**Main operations:**
- Summarisation of read abundance across the highest taxonomic groups.
- Calculation of the number of species represented within each highest taxon.
- Generation of overall and sample-level taxonomic summaries.
- Analysis of taxonomic abundance at family and species levels.
- Selection of taxonomic groups for focused analysis, including viruses, bacteria, archaea and eukaryotes.
- Identification and ranking of the most abundant taxa.
- Generation of summary tables, abundance tables and visualisations.

**Scripts:**
- `01_highest_taxon_analysis.R`
- `02_taxonomic_abundance_analysis.R`

**Main outputs:**
- Taxonomic summary and abundance tables.
- Overall and sample-level taxonomic visualisations.
- Doughnut and stacked bar charts.
- `analysis_report.csv`

### Phase 2: Taxon-Specific Read Extraction and Reference-Based Genome Assembly

This phase retrieves sequencing reads assigned to selected taxa and uses them for reference-based genome recovery.

**Main operations:**

**1. Taxon-specific read extraction**
- Extraction of reads assigned to a selected taxon and its descendant taxa using Kraken2 classification results and KrakenTools.
- Retrieval of corresponding paired-end reads from FASTQ files.
- Validation of paired-end inputs and checking of previously completed outputs.
- Generation of an extraction summary.

**2. Reference-based genome assembly**
- Mapping of extracted reads against selected reference genomes using BWA.
- Processing of mapped reads using SAMtools.
- Consensus sequence generation using iVar.
- Application of a minimum base frequency threshold of 0.5 and minimum read depth of 10 during consensus calling.
- Calculation of resolved bases and percentage of the reference genome resolved.
- Generation of consensus FASTA files and assembly summaries.

**Scripts:**
- `step_1_Read_extraction_species_with_children.sh`
- `step_2_Reference_based_assembly.sh`

**Main outputs:**
- Taxon-specific extracted FASTQ files.
- `extraction_summary.csv`
- Consensus FASTA files.
- `assembly_summary_wide.csv`

### Phase 3: Sequence Preparation and Phylogenetic Analysis

This phase prepares recovered sequences for phylogenetic analysis, constructs maximum-likelihood trees and generates metadata-annotated visualisations.

**Main operations:**
- Cleaning and standardisation of nucleotide FASTA sequences.
- Removal of duplicate sequence headers and correction of invalid nucleotide characters.
- Sequence reorientation where required.
- Multiple sequence alignment using MAFFT.
- Alignment trimming using TrimAl.
- Alignment consistency checks.
- Maximum-likelihood phylogenetic tree construction using IQ-TREE, with model selection and branch support assessment.
- Tree rerooting using a specified outgroup and Newick Utilities.
- Matching of alignment sequence identifiers to sample metadata.
- Phylogenetic tree visualisation using ggtree.
- Optional identification and extraction of closely related sequences based on phylogenetic distances.

**Scripts:**
- `01_prepare_sequences.sh`
- `02_align_and_trim.sh`
- `03_iqtree_and_reroot.sh`
- `04_extract_alignment_metadata.R`
- `05_ggtree_visualisation_edited.R`
- `close_sequences.R`

**Main outputs:**
- Cleaned and aligned FASTA files.
- Trimmed sequence alignments.
- Maximum-likelihood and consensus phylogenetic trees.
- Rooted phylogenetic trees.
- Alignment-specific metadata tables.
- Phylogenetic visualisations.
- Optional lists and FASTA files of closely related sequences.

## Repository Structure

The repository is organised into numbered directories corresponding to the four workflow steps.

```text
metagenomics-analysis-pipeline/
│
├── 01_quality_control_taxonomic_profiling/
│   ├── step_1_QC_bracken_pipeline.sh
│   ├── step_2_annotation_and_combination.py
│   └── README.md
│
├── 01b_taxonomic_abundance_analysis/
│   ├── 01_highest_taxon_analysis.R
│   ├── 02_taxonomic_abundance_analysis.R
│   └── README.md
│
├── 02_taxon_specific_read_extraction_assembly/
│   ├── step_1_Read_extraction_species_with_children.sh
│   ├── step_2_Reference_based_assembly.sh
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

## Software Requirements

The pipeline uses the following command-line tools and software:

| Software | Main application |
|---|---|
| fastp | Read quality control and filtering |
| Bowtie2 | Host-read depletion |
| Kraken2 | Taxonomic classification |
| Bracken | Taxonomic abundance estimation |
| KrakenTools | Taxon-specific read extraction |
| BWA | Reference mapping |
| SAMtools | Alignment processing |
| iVar | Consensus sequence generation |
| MAFFT | Multiple sequence alignment |
| TrimAl | Alignment trimming |
| IQ-TREE | Maximum-likelihood phylogenetic analysis |
| Newick Utilities | Phylogenetic tree rerooting |
| Python | Taxonomic annotation and data processing |
| R | Statistical summaries and visualisation |

The R scripts also require their respective R packages, including tidyverse, viridis, scales, stringr, Polychrome, ape, ggtree, ggplot2, dplyr and RColorBrewer, depending on the script.

## Getting Started

1. Clone or download the repository.
2. Install the software dependencies required for the relevant workflow step.
3. Review the README file within the directory for the step you intend to run.
4. Open each script and modify its user settings before execution.
5. Run the scripts in the appropriate workflow order, ensuring that the required input files and outputs from preceding steps are available.

### Important: Configure User Settings Before Running Scripts

**The user settings in each script must be modified to match your own working environment and input data.**

Before running any script, review its user settings section and update the relevant parameters, including:

- Project and working directory paths.
- Input FASTQ, FASTA, CSV or other required file paths.
- Output directory paths and file names, where configurable.
- Reference genome or reference sequence paths.
- Kraken2 database and NCBI taxonomy file paths, where applicable.
- Metadata file paths and sequence identifier columns.
- Target taxonomic groups or species.
- Outgroup sequence identifiers for phylogenetic tree rerooting.
- Thread numbers and other computational settings.

The exact settings will differ between scripts. Some scripts may require reference files, metadata files, database locations or specific sequence identifiers, while others may only require input and output directories.

Do not assume that the paths or reference files already specified in the scripts will work on another computer. Update them to reflect the location and organisation of your own files.

Also ensure that the output from one step is available at the location expected by the next step before proceeding.

## Reproducibility and Interpretation

The workflow is designed to support reproducible processing and analysis of metagenomic sequencing data. Software versions, reference databases, genome accessions and analysis parameters should be recorded when conducting an analysis.

Taxonomic classification and read abundance estimates indicate the assignment of sequencing reads to particular taxa. They do not, on their own, establish active infection. Similarly, genome recovery depends on read depth, sequence quality, reference availability and genome coverage.

## Citation

If you use or adapt this workflow in a research project or publication, please cite this repository using its associated GitHub or archived DOI reference, once available.
