# Phase 1b: Taxonomic Composition and Abundance Analysis

## Overview
This phase performs descriptive analysis and visualisation of taxonomic abundance estimates generated during Phase 1. The first script summarises the overall and sample-level distribution of major taxonomic groups (Highest_taxon). The second provides a more detailed, configurable analysis of taxonomic abundance at family or species level, allowing analysis of individual taxonomic groups or all groups combined.

## Scripts
- `01_highest_taxon_analysis.R`: Calculates read abundance and relative percentages for major taxonomic groups, generating summary tables, a doughnut chart and a sample-level stacked bar chart.
- `02_taxonomic_abundance_analysis.R`: Generates overall and sample-level taxonomic abundance tables, top-taxa summaries, publication-quality figures and an analysis report. Taxonomic group, taxonomic level and the number of taxa displayed can be configured.

## Requirements
R (version 4.0 or later recommended) and the following packages: `tidyverse`, `viridis`, `scales`, `stringr` and `Polychrome`. The second script automatically installs missing packages.

## Input
The consolidated species abundance dataset (`combined_species_abundance.csv`) generated in Phase 1, containing `Sample_ID`, `Highest_taxon`, `Family`, `Species` and `new_est_reads` columns.

## Execution
Update the project directory and input file paths in each script before running.

```bash
Rscript 01_highest_taxon_analysis.R

Rscript 02_taxonomic_abundance_analysis.R
```

## Outputs
Results are organised into separate `tables/` and `figures/` directories. Outputs include overall and sample-level abundance tables, top-taxa summaries, taxonomic composition figures and an `analysis_report.csv` documenting the analysis settings and summary statistics.
