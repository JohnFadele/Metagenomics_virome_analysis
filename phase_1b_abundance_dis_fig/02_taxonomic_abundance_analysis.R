###############################################################
# 02_taxonomic_abundance_analysis.R
#
# Purpose:
# Generic taxonomic abundance analysis of Bracken outputs
#
# Input:
# combined_species_abundance.csv
#
# Outputs:
#   • Summary tables
#   • Publication-quality figures
#   • Analysis report
#
###############################################################

###############################################################
## PART 1 - SETUP & USER SETTINGS
###############################################################

############################
## Load required packages
############################

required_packages <- c(
  "tidyverse",
  "viridis",
  "scales",
  "stringr",
  "Polychrome",
  "grid"
)

missing_packages <- required_packages[
  !required_packages %in% installed.packages()[, "Package"]
]

if(length(missing_packages) > 0){

  install.packages(
    missing_packages,
    dependencies = TRUE
  )

}

invisible(
  lapply(
    required_packages,
    library,
    character.only = TRUE
  )
)

############################
## User settings
############################

PROJECT_DIR <- "/home/fadele/GPADZS_3/second_phase/results"

INPUT_FILE <- file.path(
  PROJECT_DIR,
  "gpazds_3_combined_species_abundance.csv"
)

#################################################
## Analysis settings
#################################################

# Options:
# "All"
# "Bacteria"
# "Viruses"
# "Archaea"
# "Eukaryota"
# or any value present in Highest_taxon

TARGET_TAXON <- "Viruses"

#################################################

# Options:
# "Family"
# "Species"

TARGET_LEVEL <- "Family"

#################################################

TOP_OVERALL <- 20

TOP_SAMPLE <- 10

############################
## Figure settings
############################

DPI <- 800

WIDTH_OVERALL <- 12
HEIGHT_OVERALL <- 8

WIDTH_SAMPLE <- 32
HEIGHT_SAMPLE <- 17

BASE_FONT_SIZE <- 13

LEGEND_ROWS <- 9

############################
## Output folders
############################

OUTPUT_DIR <- file.path(

  PROJECT_DIR,

  "descriptive_analysis",

  paste0(

    tolower(TARGET_TAXON),

    "_",

    tolower(TARGET_LEVEL)

  )

)

TABLE_DIR <- file.path(
  OUTPUT_DIR,
  "tables"
)

FIGURE_DIR <- file.path(
  OUTPUT_DIR,
  "figures"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

############################
## Validate input file
############################

if(!file.exists(INPUT_FILE)){

  stop(
    paste(
      "\nInput file not found:\n",
      INPUT_FILE
    )
  )

}

############################
## Start analysis
############################

cat("\n")
cat("=========================================================\n")
cat("           TAXONOMIC ABUNDANCE ANALYSIS\n")
cat("=========================================================\n")
cat("Input file      :", INPUT_FILE, "\n")
cat("Target taxon    :", TARGET_TAXON, "\n")
cat("Taxonomic level :", TARGET_LEVEL, "\n")
cat("Output folder   :", OUTPUT_DIR, "\n")
cat("=========================================================\n\n")

###############################################################
## PART 2 - IMPORT, VALIDATE & PREPARE DATA
###############################################################

cat("Reading input data...\n")

############################
## Read input file
############################

dat <- read_csv(
  INPUT_FILE,
  show_col_types = FALSE
)

cat("Input file successfully loaded.\n\n")

############################
## Validate required columns
############################

required_columns <- c(

  "Sample_ID",
  "Highest_taxon",
  "Family",
  "Species",
  "new_est_reads"

)

missing_columns <- setdiff(
  required_columns,
  names(dat)
)

if(length(missing_columns) > 0){

  stop(

    paste(

      "The following required columns are missing:\n\n",

      paste(missing_columns, collapse = ", ")

    )

  )

}

############################
## Convert read counts
############################

dat <- dat %>%

  mutate(

    new_est_reads = as.numeric(new_est_reads)

  )

############################
## Remove invalid records
############################

dat <- dat %>%

  filter(

    !is.na(new_est_reads),

    new_est_reads > 0

  )

############################
## Display available taxa
############################

available_taxa <- dat %>%

  distinct(Highest_taxon) %>%

  arrange(Highest_taxon) %>%

  pull(Highest_taxon)

cat("Available Highest_taxon values:\n")

cat(
  paste(
    available_taxa,
    collapse = "\n"
  )
)

cat("\n\n")

############################
## Filter selected taxon
############################

if(TARGET_TAXON != "All"){

  if(!(TARGET_TAXON %in% available_taxa)){

    stop(

      paste(

        "\nRequested Highest_taxon not found.\n\n",

        "Available values are:\n",

        paste(
          available_taxa,
          collapse = "\n"
        )

      )

    )

  }

  dat <- dat %>%

    filter(

      Highest_taxon == TARGET_TAXON

    )

}

############################
## Select taxonomic level
############################

if(TARGET_LEVEL == "Family"){

  dat <- dat %>%

    mutate(

      Taxon = Family

    )

}else if(TARGET_LEVEL == "Species"){

  dat <- dat %>%

    mutate(

      Taxon = Species

    )

}else{

  stop(

    "TARGET_LEVEL must be either 'Family' or 'Species'."

  )

}

############################
## Remove missing taxa
############################

dat <- dat %>%

  filter(

    !is.na(Taxon),

    Taxon != "",

    Taxon != "Unknown",

    Taxon != "unclassified",

    Taxon != "Unclassified"

  )

############################
## Standardise sample order
############################

sample_order <- dat %>%

  distinct(Sample_ID) %>%

  arrange(Sample_ID) %>%

  pull(Sample_ID)

dat <- dat %>%

  mutate(

    Sample_ID = factor(

      Sample_ID,

      levels = sample_order

    )

  )

############################
## Basic statistics
############################

TOTAL_READS <- sum(dat$new_est_reads)

TOTAL_SAMPLES <- n_distinct(dat$Sample_ID)

TOTAL_TAXA <- n_distinct(dat$Taxon)

############################
## Check dataset
############################

if(nrow(dat) == 0){

  stop(

    "No records remain after filtering."

  )

}

###############################################################
## Analysis summary
###############################################################

cat("=========================================\n")

cat("Analysis dataset prepared successfully.\n\n")

cat("Target taxon    :", TARGET_TAXON, "\n")

cat("Taxonomic level :", TARGET_LEVEL, "\n")

cat("Samples         :", TOTAL_SAMPLES, "\n")

cat("Unique taxa     :", TOTAL_TAXA, "\n")

cat("Total reads     :", comma(TOTAL_READS), "\n")

cat("=========================================\n\n")

###############################################################
## PART 3 - GENERATE SUMMARY TABLES
###############################################################

cat("Generating summary tables...\n\n")

###############################################################
## Overall abundance
###############################################################

overall_abundance <-

  dat %>%

  group_by(Taxon) %>%

  summarise(

    Reads = sum(new_est_reads),

    Number_of_Samples = n_distinct(Sample_ID),

    .groups = "drop"

  ) %>%

  mutate(

    Percentage = (Reads / sum(Reads)) * 100

  ) %>%

  arrange(

    desc(Reads)

  ) %>%

  mutate(

    Rank = row_number()

  )

###############################################################
## Sample abundance
###############################################################

sample_abundance <-

  dat %>%

  group_by(

    Sample_ID,

    Taxon

  ) %>%

  summarise(

    Reads = sum(new_est_reads),

    .groups = "drop"

  ) %>%

  group_by(

    Sample_ID

  ) %>%

  mutate(

    Percentage = (Reads / sum(Reads)) * 100

  ) %>%

  ungroup()

###############################################################
## Overall Top N
###############################################################

top_abundance <-

  overall_abundance %>%

  slice_head(

    n = TOP_OVERALL

  )

###############################################################
## Summary statistics
###############################################################

if(nrow(overall_abundance) > 0){

  summary_table <-

    tibble(

      Target_Taxon = TARGET_TAXON,

      Taxonomic_Level = TARGET_LEVEL,

      Samples = TOTAL_SAMPLES,

      Unique_Taxa = TOTAL_TAXA,

      Total_Reads = TOTAL_READS,

      Top_Taxon = overall_abundance$Taxon[1],

      Top_Taxon_Reads = overall_abundance$Reads[1],

      Top_Taxon_Percentage = round(

        overall_abundance$Percentage[1],

        2

      )

    )

}else{

  summary_table <-

    tibble(

      Target_Taxon = TARGET_TAXON,

      Taxonomic_Level = TARGET_LEVEL,

      Samples = TOTAL_SAMPLES,

      Unique_Taxa = TOTAL_TAXA,

      Total_Reads = TOTAL_READS,

      Top_Taxon = NA_character_,

      Top_Taxon_Reads = NA_real_,

      Top_Taxon_Percentage = NA_real_

    )

}

###############################################################
## Export tables
###############################################################

write_csv(

  overall_abundance,

  file.path(

    TABLE_DIR,

    paste0(

      "overall_",

      tolower(TARGET_TAXON),

      "_",

      tolower(TARGET_LEVEL),

      ".csv"

    )

  )

)

write_csv(

  sample_abundance,

  file.path(

    TABLE_DIR,

    paste0(

      "sample_",

      tolower(TARGET_TAXON),

      "_",

      tolower(TARGET_LEVEL),

      ".csv"

    )

  )

)

write_csv(

  top_abundance,

  file.path(

    TABLE_DIR,

    paste0(

      "top",

      TOP_OVERALL,

      "_",

      tolower(TARGET_TAXON),

      "_",

      tolower(TARGET_LEVEL),

      ".csv"

    )

  )

)

write_csv(

  summary_table,

  file.path(

    TABLE_DIR,

    paste0(

      "summary_",

      tolower(TARGET_TAXON),

      "_",

      tolower(TARGET_LEVEL),

      ".csv"

    )

  )

)

###############################################################
## Console report
###############################################################

cat("Summary tables generated successfully.\n\n")

cat("Unique taxa :", nrow(overall_abundance), "\n")

cat("Top taxa    :", nrow(top_abundance), "\n\n")


###############################################################
## PART 4 - PREPARE PLOTTING DATA
###############################################################

cat("Preparing plotting datasets...\n\n")

###############################################################
## Publication theme
###############################################################

publication_theme <-

  theme_bw(base_size = BASE_FONT_SIZE) +

  theme(

    panel.grid.major = element_blank(),

    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      colour = "black",
      linewidth = 0.7
    ),

    axis.title = element_text(
      face = "bold"
    ),

    axis.text = element_text(
      colour = "black"
    ),

    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),

    legend.position = "bottom",

    legend.title = element_blank(),

    legend.text = element_text(
      size = 10
    ),

    legend.key.width = unit(
      0.6,
      "cm"
    ),

    legend.key.height = unit(
      0.4,
      "cm"
    ),

    legend.spacing.x = unit(
      0.2,
      "cm"
    )

  )

###############################################################
## Function for saving publication figures
###############################################################

save_publication_plot <- function(
  plot,
  filename,
  width,
  height
){

  ggsave(

    filename = filename,

    plot = plot,

    width = width,

    height = height,

    dpi = DPI,

    units = "in",

    limitsize = FALSE

  )

}

###############################################################
## Figure 1 dataset
###############################################################

overall_plot <-

  top_abundance %>%

  mutate(

    Taxon = factor(

      Taxon,

      levels = rev(Taxon)

    )

  )

###############################################################
## Figure 2 dataset
###############################################################

sample_plot <-

  sample_abundance %>%

  group_by(
    Sample_ID
  ) %>%

  arrange(
    desc(Percentage),
    .by_group = TRUE
  ) %>%

  mutate(

    Rank = row_number(),

    Plot_Taxon = if_else(

      Rank <= TOP_SAMPLE,

      Taxon,

      "Others"

    )

  ) %>%

  group_by(
    Sample_ID,
    Plot_Taxon
  ) %>%

  summarise(

    Percentage = sum(Percentage),

    .groups = "drop"

  ) %>%

  ungroup()

###############################################################
## Preserve sample order
###############################################################

sample_plot <-

  sample_plot %>%

  mutate(

    Sample_ID = factor(

      Sample_ID,

      levels = levels(dat$Sample_ID)

    )

  )

###############################################################
## Colour palette
###############################################################

plot_taxa <-

  sample_plot %>%

  filter(
    Plot_Taxon != "Others"
  ) %>%

  group_by(
    Plot_Taxon
  ) %>%

  summarise(

    Total = sum(Percentage),

    .groups = "drop"

  ) %>%

  arrange(
    desc(Total)
  ) %>%

  pull(
    Plot_Taxon
  )

###############################################################
## Generate colours automatically
###############################################################

plot_colours <-

  setNames(

    Polychrome::createPalette(

      N = length(plot_taxa),

      seedcolors = Polychrome::palette36.colors(8)

    ),

    plot_taxa

  )

###############################################################
## Others
###############################################################

plot_colours <-

  c(

    plot_colours,

    Others = "grey80"

  )

###############################################################
## Legend order
###############################################################

sample_plot <-

  sample_plot %>%

  mutate(

    Plot_Taxon = factor(

      Plot_Taxon,

      levels = c(
        plot_taxa,
        "Others"
      )

    )

  )

###############################################################
## Console report
###############################################################

cat("Plotting datasets prepared.\n\n")

cat("Overall taxa : ", nrow(overall_plot), "\n")

cat("Legend taxa  : ", length(plot_taxa), "\n")

cat("Colours      : ", length(plot_colours), "\n\n")

###############################################################
## PART 5 - GENERATE FIGURES
###############################################################

cat("Generating figures...\n\n")

###############################################################
## Figure titles
###############################################################

if (TARGET_TAXON == "All") {

  OVERALL_TITLE <- paste(
    "Top", TOP_OVERALL, TARGET_LEVEL
  )

  SAMPLE_TITLE <- paste(
    "Top", TOP_SAMPLE, TARGET_LEVEL,
    "within each sample"
  )

} else {

  taxon_name <- switch(

    TARGET_TAXON,

    "Bacteria" = "Bacterial",
    "Viruses" = "Viral",
    "Archaea" = "Archaeal",
    "Eukaryota" = "Eukaryotic",

    TARGET_TAXON

  )

  OVERALL_TITLE <- paste(
    "Top",
    TOP_OVERALL,
    taxon_name,
    TARGET_LEVEL
  )

  SAMPLE_TITLE <- paste(
    "Top",
    TOP_SAMPLE,
    taxon_name,
    TARGET_LEVEL,
    "within each sample"
  )

}

###############################################################
## Figure 1
## Overall abundance
###############################################################

figure1 <-

  ggplot(

    overall_plot,

    aes(

      x = Percentage,

      y = Taxon,

      fill = Taxon

    )

  ) +

  geom_col(

    colour = "black",

    linewidth = 0.2

  ) +

  scale_fill_manual(

    values = plot_colours,

    drop = FALSE

  ) +

  scale_x_continuous(

    labels = percent_format(scale = 1),

    expand = expansion(mult = c(0, 0.02))

  ) +

  labs(

    title = OVERALL_TITLE,

    x = "Relative abundance (%)",

    y = NULL

  ) +

  publication_theme +

  theme(

    legend.position = "none"

  )

###############################################################
## Save Figure 1
###############################################################

save_publication_plot(

  plot = figure1,

  filename = file.path(

    FIGURE_DIR,

    paste0(

      "overall_",

      tolower(TARGET_TAXON),

      "_",

      tolower(TARGET_LEVEL),

      ".png"

    )

  ),

  width = WIDTH_OVERALL,

  height = HEIGHT_OVERALL

)

###############################################################
## Figure 2
###############################################################

figure2 <-

  ggplot(

    sample_plot,

    aes(

      x = Sample_ID,

      y = Percentage,

      fill = Plot_Taxon

    )

  ) +

  geom_col(

    colour = "black",

    linewidth = 0.15

  ) +

  scale_fill_manual(

    values = plot_colours,

    drop = FALSE

  ) +

  scale_y_continuous(

    labels = percent_format(scale = 1),

    expand = expansion(mult = c(0, 0))

  ) +

  labs(

    title = SAMPLE_TITLE,

    x = "Sample",

    y = "Relative abundance (%)"

  ) +

  publication_theme +

  guides(

    fill = guide_legend(

      nrow = LEGEND_ROWS,

      byrow = TRUE

    )

  ) +

  theme(

    axis.text.x = element_text(

      angle = 90,

      hjust = 1,

      vjust = 0.5

    )

  )

###############################################################
## Save Figure 2
###############################################################

save_publication_plot(

  plot = figure2,

  filename = file.path(

    FIGURE_DIR,

    paste0(

      "sample_",

      tolower(TARGET_TAXON),

      "_",

      tolower(TARGET_LEVEL),

      ".png"

    )

  ),

  width = WIDTH_SAMPLE,

  height = HEIGHT_SAMPLE

)

###############################################################
## Report
###############################################################

cat("Figures generated successfully.\n\n")

cat("Figure 1 saved.\n")

cat("Figure 2 saved.\n\n")

###############################################################
## PART 6 - ANALYSIS REPORT
###############################################################

cat("Generating analysis report...\n\n")

###############################################################
## Create run report
###############################################################

run_report <-

  tibble(

    Analysis_Date = Sys.Date(),

    Analysis_Time = format(Sys.time(), "%H:%M:%S"),

    Input_File = INPUT_FILE,

    Target_Taxon = TARGET_TAXON,

    Taxonomic_Level = TARGET_LEVEL,

    Samples = TOTAL_SAMPLES,

    Unique_Taxa = TOTAL_TAXA,

    Total_Reads = TOTAL_READS,

    Top_Overall_Taxa = nrow(top_abundance),

    Top_Sample_Taxa = TOP_SAMPLE,

    Output_Directory = OUTPUT_DIR

  )

write_csv(

  run_report,

  file.path(

    OUTPUT_DIR,

    "analysis_report.csv"

  )

)

###############################################################
## Completion message
###############################################################

cat("=====================================================\n")

cat("Analysis completed successfully.\n\n")

cat("Input file:\n")

cat(INPUT_FILE, "\n\n")

cat("Output directory:\n")

cat(OUTPUT_DIR, "\n\n")

cat("Summary\n")

cat("---------------------------------------------\n")

cat("Samples analysed :", TOTAL_SAMPLES, "\n")

cat("Unique taxa      :", TOTAL_TAXA, "\n")

cat("Total reads      :", comma(TOTAL_READS), "\n")

cat("---------------------------------------------\n\n")

cat("Tables\n")

cat("---------------------------------------------\n")

cat("overall_", tolower(TARGET_TAXON), "_",
    tolower(TARGET_LEVEL), ".csv\n", sep="")

cat("sample_", tolower(TARGET_TAXON), "_",
    tolower(TARGET_LEVEL), ".csv\n", sep="")

cat("top", TOP_OVERALL, "_",
    tolower(TARGET_TAXON), "_",
    tolower(TARGET_LEVEL), ".csv\n", sep="")

cat("summary_", tolower(TARGET_TAXON), "_",
    tolower(TARGET_LEVEL), ".csv\n", sep="")

cat("---------------------------------------------\n\n")

cat("Figures\n")

cat("---------------------------------------------\n")

cat("overall_", tolower(TARGET_TAXON), "_",
    tolower(TARGET_LEVEL), ".png\n", sep="")

cat("sample_", tolower(TARGET_TAXON), "_",
    tolower(TARGET_LEVEL), ".png\n", sep="")

cat("---------------------------------------------\n\n")

cat("Run report\n")

cat("---------------------------------------------\n")

cat("analysis_report.csv\n\n")

cat("=====================================================\n")





