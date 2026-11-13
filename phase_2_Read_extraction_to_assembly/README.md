# Phase 2: Taxon-Specific Read Extraction and Reference-Based Genome Assembly

## Overview
This phase extracts paired-end sequencing reads assigned to a specified taxonomic identifier and all its descendant taxa using Kraken2 classification outputs and KrakenTools. The extracted reads are subsequently mapped against reference genomes using BWA-MEM, with consensus sequences generated using iVar. Assembly completeness is assessed by calculating the number and percentage of resolved reference bases.

## Scripts
- `step_1_Read_extraction_species_with_children.sh`: Extracts reads assigned to a target taxon and all its descendants from paired-end FASTQ files. Includes sample-pair validation, extraction status tracking and support for resuming completed samples.
- `step_2_Reference_based_assembly.sh`: Performs reference-based mapping and consensus generation across multiple samples and reference genomes. Produces consensus FASTA files, calculates resolved bases and percentage reference coverage, and generates a consolidated assembly summary.

## Requirements
- **Software:** Python 3, KrakenTools, BWA, SAMtools and iVar.
- **Input resources:** Kraken2 classification outputs and reports, paired-end FASTQ files, and reference genome FASTA files.

## Input
Paired-end FASTQ files (`.fastq` or `.fastq.gz`) and corresponding Kraken2 classification outputs and reports. Reference genomes must be provided in FASTA format (`.fa`, `.fasta` or `.fna`).

## Execution
Update the directory paths, target taxonomic identifier, target name and computational settings in both scripts before execution.

Run sequentially:

```bash
bash step_1_Read_extraction_species_with_children.sh

bash step_2_Reference_based_assembly.sh
```

## Outputs
Read extraction produces taxon-specific paired-end FASTQ files and an `extraction_summary.csv`. Reference-based assembly generates sample-specific consensus FASTA files, a central collection of consensus sequences, processing logs and an `assembly_summary_wide.csv` containing resolved bases, reference lengths and percentage reference resolution for each sample-reference combination.
