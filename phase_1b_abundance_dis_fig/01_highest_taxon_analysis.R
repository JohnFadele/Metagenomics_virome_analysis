###############################################################
## 01_highest_taxon_analysis.R
##
## Calculates and visualises the relative abundance of the
## major taxonomic groups (Highest_taxon)
##
## Input:
## combined_species_abundance.csv
##
## Outputs:
## Tables:
##   overall_highest_taxon.csv
##   overall_highest_taxon_summary.csv
##   sample_highest_taxon.csv
##
## Figures:
##   overall_highest_taxon_doughnut.png
##   sample_highest_taxon_stacked_bar.png
###############################################################

#############################
## Libraries
#############################

library(tidyverse)
library(viridis)
library(scales)

#############################
## User settings
#############################

PROJECT_DIR <- "/home/fadele/GPADZS_3/second_phase/results"

INPUT_FILE <- file.path(
  PROJECT_DIR,
  "gpazds_3_combined_species_abundance.csv"
)

OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "descriptive_analysis"
)

TABLE_DIR <- file.path(OUTPUT_DIR,"tables")
FIGURE_DIR <- file.path(OUTPUT_DIR,"figures")

dir.create(TABLE_DIR,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(FIGURE_DIR,
           recursive = TRUE,
           showWarnings = FALSE)

#############################
## Read data
#############################

dat <- read.csv(
  INPUT_FILE,
  stringsAsFactors = FALSE
)

#############################
## Data cleaning
#############################

dat <- dat %>%
  filter(Highest_taxon != "") %>%
  mutate(
    new_est_reads = as.numeric(new_est_reads)
  )

#############################
## Overall composition
#############################

overall <- dat %>%
  group_by(Highest_taxon) %>%
  summarise(
    Reads = sum(new_est_reads),
    Number_of_species = n_distinct(Species),
    .groups = "drop"
  ) %>%
  mutate(
    Percentage = Reads/sum(Reads)*100
  ) %>%
  arrange(desc(Percentage))

write.csv(
  overall,
  file.path(
    TABLE_DIR,
    "overall_highest_taxon.csv"
  ),
  row.names = FALSE
)

#############################
## Overall summary
#############################

overall_summary <- data.frame(

  Total_reads =
    sum(overall$Reads),

  Total_highest_taxa =
    nrow(overall),

  Total_species =
    n_distinct(dat$Species)

)

write.csv(

  overall_summary,

  file.path(
    TABLE_DIR,
    "overall_highest_taxon_summary.csv"
  ),

  row.names = FALSE

)

#############################
## Sample composition
#############################

sample_summary <- dat %>%

  group_by(
    Sample_ID,
    Highest_taxon
  ) %>%

  summarise(

    Reads =
      sum(new_est_reads),

    .groups="drop"

  ) %>%

  group_by(
    Sample_ID
  ) %>%

  mutate(

    Percentage =
      Reads/sum(Reads)*100

  ) %>%

  ungroup()

sample_summary$Sample_ID <-

  factor(

    sample_summary$Sample_ID,

    levels=
      sort(unique(sample_summary$Sample_ID))

  )

write.csv(

  sample_summary,

  file.path(
    TABLE_DIR,
    "sample_highest_taxon.csv"
  ),

  row.names=FALSE

)

###########################################################
## Figure 1
## Doughnut
###########################################################

total_reads <-

  comma(sum(overall$Reads))

p1 <-

ggplot(

  overall,

  aes(
    x=2,
    y=Percentage,
    fill=Highest_taxon
  )

)+

geom_col(

  colour="white",
  linewidth=.7

)+

coord_polar(theta="y")+

xlim(.5,2.5)+

geom_text(

  aes(

    label=paste0(

      round(Percentage,1),

      "%"

    )

  ),

  position=
    position_stack(vjust=.5),

  size=5

)+

annotate(

  "text",

  x=0,

  y=0,

  label=paste(

    "Total reads\n",

    total_reads

  ),

  fontface="bold",

  size=6

)+

scale_fill_viridis_d(option="D")+

theme_void(base_size=16)+

theme(

  legend.position="right",

  plot.title=
    element_text(
      hjust=.5,
      face="bold"
    )

)+

labs(

  title=
    "Overall Highest Taxonomic Composition",

  fill="Highest taxon"

)

ggsave(

  filename=file.path(

    FIGURE_DIR,

    "overall_highest_taxon_doughnut.png"

  ),

  plot=p1,

  dpi=1000,

  width=8,

  height=8

)

###########################################################
## Figure 2
## 100% stacked bar
###########################################################

p2 <-

ggplot(

  sample_summary,

  aes(

    Sample_ID,

    Percentage,

    fill=Highest_taxon

  )

)+

geom_bar(

  stat="identity",

  width=.8

)+

scale_fill_viridis_d(option="D")+

scale_y_continuous(

  labels=function(x) paste0(x,"%"),

  expand=c(0,0)

)+

labs(

  x="Sample",

  y="Relative abundance (%)",

  fill="Highest taxon",

  title="Highest Taxonomic Composition per Sample"

)+

theme_bw(base_size=15)+

theme(

  plot.title=
    element_text(
      face="bold",
      hjust=.5
    ),

  axis.text.x=
    element_text(
      angle=90,
      hjust=1
    ),

  panel.grid.minor=
    element_blank(),

  panel.grid.major.x=
    element_blank()

)

ggsave(

  filename=file.path(

    FIGURE_DIR,

    "sample_highest_taxon_stacked_bar.png"

  ),

  plot=p2,

  dpi=1000,

  width=25,

  height=7

)

###########################################################

cat("\nAnalysis completed successfully.\n")
