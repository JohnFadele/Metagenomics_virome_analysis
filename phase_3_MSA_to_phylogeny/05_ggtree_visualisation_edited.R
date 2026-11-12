# ============================================================
# PPR PHYLOGENETIC TREE VISUALISATION
# ============================================================


# ============================================================
# USER SETTINGS
# ============================================================

PROJECT_DIR <- "/home/fadele/ARS_NACA/ppr_phylo"

TREE_FILE <- file.path(
  PROJECT_DIR,
  "03_IQTREE",
  "ppr_pub_study_og_rooted.treefile"
)

METADATA_FILE <- file.path(
  PROJECT_DIR,
  "04_metadata",
  "ppr_pub_study_og_alignment_metadata.csv"
)

OUTPUT_FILE <- file.path(
  PROJECT_DIR,
  "05_visualisation",
  "PPR_WG_phylogeny_genotype_country.png"
)

DROP_TIPS <- c(
  "NC_006296.2"
)

SEQUENCE_ID_COLUMN <- "sequence_id"
# Metadata column used to colour the tree tips.
# Change ONLY this value when you want to colour tips by a different
# metadata variable, e.g. "Country", "Region", "Continent", "Host", etc.
TIP_COLOUR_COLUMN <- "Region"

# Number of rows used for the tip-colour legend.
TIP_LEGEND_NROW <- 1

COUNTRY_COLUMN <- "Country"
REGION_COLUMN <- "Region"
GENOTYPE_COLUMN <- "Genotype"
STUDY_SAMPLE_COLUMN <- "study_sample"
STUDY_SAMPLE_VALUE <- "yes"


# ============================================================
# GENOTYPE VISUALISATION MODE
# ============================================================

# "all"        = retain all sequences
# "known_only" = retain only sequences with known genotype

GENOTYPE_MODE <- "all"


# ============================================================
# GENOTYPE HIGHLIGHT SETTINGS
# ============================================================

MIN_KNOWN_GENOTYPE_TIPS <- 1

HIGHLIGHT_ALPHA <- 0.28

HIGHLIGHT_EXTEND_FRACTION <- 1.03


# ============================================================
# TREE APPEARANCE
# ============================================================

TREE_LINEWIDTH <- 0.45

TREE_COLOUR <- "#263238"

TIP_SIZE <- 3

# Study samples are shown as a larger filled point over the normal
# country-coloured tip. Shape 21 allows fill and border colour to be
# controlled independently.
STUDY_POINT_SHAPE <- 21
STUDY_POINT_SIZE <- 4.2
STUDY_POINT_STROKE <- 1.25


# ============================================================
# OUTPUT SETTINGS
# ============================================================

OUTPUT_WIDTH <- 12

OUTPUT_HEIGHT <- 12

OUTPUT_DPI <- 1000


# ============================================================
# CHECK GENOTYPE MODE
# ============================================================

if (!GENOTYPE_MODE %in% c("all", "known_only")) {

  stop(
    "GENOTYPE_MODE must be either 'all' or 'known_only'."
  )

}


# ============================================================
# LOAD PACKAGES
# ============================================================

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(RColorBrewer)


# ============================================================
# CREATE OUTPUT DIRECTORY
# ============================================================

OUTPUT_DIR <- dirname(OUTPUT_FILE)

if (!dir.exists(OUTPUT_DIR)) {

  dir.create(
    OUTPUT_DIR,
    recursive = TRUE
  )

}


# ============================================================
# CHECK INPUT FILES
# ============================================================

if (!file.exists(TREE_FILE)) {

  stop(
    paste(
      "Tree file not found:",
      TREE_FILE
    )
  )

}

if (!file.exists(METADATA_FILE)) {

  stop(
    paste(
      "Metadata file not found:",
      METADATA_FILE
    )
  )

}


# ============================================================
# ANALYSIS NAME
# ============================================================

analysis_name <- basename(TREE_FILE)

analysis_name <- sub(
  "_rooted\\.treefile$",
  "",
  analysis_name
)


# ============================================================
# READ PHYLOGENETIC TREE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("READING PHYLOGENETIC TREE\n")
cat("============================================================\n")

tree <- read.tree(TREE_FILE)

tree$tip.label <- trimws(
  tree$tip.label
)

cat(
  "Initial number of tips:",
  length(tree$tip.label),
  "\n"
)


# ============================================================
# REMOVE SPECIFIED TIPS
# ============================================================

if (length(DROP_TIPS) > 0) {

  tips_present <- DROP_TIPS[
    DROP_TIPS %in% tree$tip.label
  ]

  tips_missing <- DROP_TIPS[
    !DROP_TIPS %in% tree$tip.label
  ]

  if (length(tips_present) > 0) {

    tree <- drop.tip(
      tree,
      tips_present
    )

    cat(
      "Tips removed:",
      length(tips_present),
      "\n"
    )

    cat(
      "Removed:",
      paste(
        tips_present,
        collapse = ", "
      ),
      "\n"
    )

  }

  if (length(tips_missing) > 0) {

    cat(
      "Requested tips not found:",
      paste(
        tips_missing,
        collapse = ", "
      ),
      "\n"
    )

  }

}

cat(
  "Tips after explicit removal:",
  length(tree$tip.label),
  "\n"
)


# ============================================================
# READ METADATA
# ============================================================

cat("\n")
cat("============================================================\n")
cat("READING METADATA\n")
cat("============================================================\n")

metadata <- read.csv(
  METADATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  SEQUENCE_ID_COLUMN,
  COUNTRY_COLUMN,
  REGION_COLUMN,
  TIP_COLOUR_COLUMN,
  GENOTYPE_COLUMN,
  STUDY_SAMPLE_COLUMN
)

missing_columns <- required_columns[
  !required_columns %in% colnames(metadata)
]

if (length(missing_columns) > 0) {

  stop(
    paste(
      "Missing metadata columns:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )

}


# ============================================================
# CLEAN METADATA
# ============================================================

metadata[[SEQUENCE_ID_COLUMN]] <- trimws(
  as.character(
    metadata[[SEQUENCE_ID_COLUMN]]
  )
)

metadata[[COUNTRY_COLUMN]] <- trimws(
  as.character(
    metadata[[COUNTRY_COLUMN]]
  )
)

metadata[[REGION_COLUMN]] <- trimws(
  as.character(
    metadata[[REGION_COLUMN]]
  )
)

metadata[[GENOTYPE_COLUMN]] <- trimws(
  as.character(
    metadata[[GENOTYPE_COLUMN]]
  )
)

metadata[[STUDY_SAMPLE_COLUMN]] <- trimws(
  as.character(
    metadata[[STUDY_SAMPLE_COLUMN]]
  )
)


# Blank genotype = unknown genotype

metadata[[GENOTYPE_COLUMN]][
  metadata[[GENOTYPE_COLUMN]] == ""
] <- NA_character_


cat(
  "Metadata records:",
  nrow(metadata),
  "\n"
)

cat(
  "Known genotype records:",
  sum(
    !is.na(
      metadata[[GENOTYPE_COLUMN]]
    )
  ),
  "\n"
)


# ============================================================
# APPLY GENOTYPE MODE
# ============================================================

if (GENOTYPE_MODE == "known_only") {

  cat("\n")
  cat("============================================================\n")
  cat("GENOTYPE MODE: KNOWN ONLY\n")
  cat("============================================================\n")

  known_metadata <- metadata[
    !is.na(
      metadata[[GENOTYPE_COLUMN]]
    ) &
      metadata[[GENOTYPE_COLUMN]] != "",
    ,
    drop = FALSE
  ]

  # IMPORTANT:
  # Extract the sequence_id column correctly.

  known_ids <- known_metadata[[SEQUENCE_ID_COLUMN]]

  known_ids <- trimws(
    as.character(
      known_ids
    )
  )

  original_tree_tips <- tree$tip.label

  tips_to_remove <- original_tree_tips[
    !original_tree_tips %in% known_ids
  ]

  if (length(tips_to_remove) > 0) {

    tree <- drop.tip(
      tree,
      tips_to_remove
    )

  }

  metadata <- known_metadata

  cat(
    "Sequences retained:",
    length(tree$tip.label),
    "\n"
  )

  cat(
    "Sequences removed because genotype was unknown:",
    length(tips_to_remove),
    "\n"
  )

} else {

  cat("\n")
  cat("============================================================\n")
  cat("GENOTYPE MODE: ALL SEQUENCES\n")
  cat("============================================================\n")

  cat(
    "All remaining sequences will be retained.\n"
  )

}


# ============================================================
# MATCH TREE SEQUENCES TO METADATA
# ============================================================

tree_labels <- trimws(
  tree$tip.label
)

metadata_ids <- trimws(
  metadata[[SEQUENCE_ID_COLUMN]]
)

matched_indices <- match(
  tree_labels,
  metadata_ids
)

matched_count <- sum(
  !is.na(matched_indices)
)

unmatched_count <- sum(
  is.na(matched_indices)
)

cat("\n")

cat(
  "Tree sequences with matching metadata:",
  matched_count,
  "\n"
)

cat(
  "Tree sequences without matching metadata:",
  unmatched_count,
  "\n"
)


# ============================================================
# IDENTIFY ROOT NODE
# ============================================================

root_node <- setdiff(
  tree$edge[, 1],
  tree$edge[, 2]
)

if (length(root_node) != 1) {

  stop(
    "Could not identify a unique root node."
  )

}

root_node <- root_node[1]


# ============================================================
# FUNCTION: GET DESCENDANT TIPS
# ============================================================

get_descendant_tips <- function(
  tree,
  node
) {

  descendants <- node

  repeat {

    child_nodes <- tree$edge[
      tree$edge[, 1] %in% descendants,
      2
    ]

    new_nodes <- child_nodes[
      !child_nodes %in% descendants
    ]

    if (length(new_nodes) == 0) {

      break

    }

    descendants <- c(
      descendants,
      new_nodes
    )

  }

  tip_nodes <- descendants[
    descendants <= Ntip(tree)
  ]

  return(
    tree$tip.label[
      tip_nodes
    ]
  )

}


# ============================================================
# SCAN INTERNAL NODES FOR GENOTYPE CLADES
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SCANNING FOR GENOTYPE-CONSISTENT CLADES\n")
cat("============================================================\n")

candidate_clades <- data.frame(
  node = integer(),
  genotype = character(),
  total_tips = integer(),
  known_tips = integer(),
  stringsAsFactors = FALSE
)

internal_nodes <- unique(
  tree$edge[, 1]
)

internal_nodes <- internal_nodes[
  internal_nodes > Ntip(tree)
]


# ============================================================
# EXAMINE EACH INTERNAL NODE
# ============================================================

for (node in internal_nodes) {

  if (node == root_node) {

    next

  }

  descendant_labels <- get_descendant_tips(
    tree,
    node
  )

  if (length(descendant_labels) == 0) {

    next

  }

  descendant_metadata <- metadata[
    match(
      descendant_labels,
      metadata[[SEQUENCE_ID_COLUMN]]
    ),
    ,
    drop = FALSE
  ]

  # THIS IS THE WORKING GENOTYPE EXTRACTION.
  # DO NOT CHANGE THIS.

  genotype_values <- descendant_metadata[[GENOTYPE_COLUMN]]

  genotype_values <- genotype_values[
    !is.na(genotype_values) &
      genotype_values != ""
  ]

  known_count <- length(
    genotype_values
  )

  if (
    known_count <
    MIN_KNOWN_GENOTYPE_TIPS
  ) {

    next

  }

  unique_genotypes <- unique(
    genotype_values
  )

  # Only highlight a clade where all known
  # genotype records agree.

  if (length(unique_genotypes) != 1) {

    next

  }

  candidate_clades <- rbind(
    candidate_clades,
    data.frame(
      node = node,
      genotype = unique_genotypes[1],
      total_tips = length(
        descendant_labels
      ),
      known_tips = known_count,
      stringsAsFactors = FALSE
    )
  )

}


# ============================================================
# REMOVE NESTED / REDUNDANT HIGHLIGHTS
# ============================================================

if (nrow(candidate_clades) > 0) {

  candidate_clades <- candidate_clades[
    order(
      -candidate_clades$total_tips,
      -candidate_clades$known_tips
    ),
    ,
    drop = FALSE
  ]

  selected_clades <- data.frame(
    node = integer(),
    genotype = character(),
    total_tips = integer(),
    known_tips = integer(),
    stringsAsFactors = FALSE
  )

  selected_nodes <- integer()

  for (
    i in seq_len(
      nrow(candidate_clades)
    )
  ) {

    current_node <- candidate_clades$node[i]

    current_descendants <- get_descendant_tips(
      tree,
      current_node
    )

    already_covered <- FALSE

    if (length(selected_nodes) > 0) {

      for (
        selected_node in selected_nodes
      ) {

        selected_descendants <- get_descendant_tips(
          tree,
          selected_node
        )

        if (
          all(
            current_descendants %in%
              selected_descendants
          )
        ) {

          already_covered <- TRUE

          break

        }

      }

    }

    if (!already_covered) {

      selected_clades <- rbind(
        selected_clades,
        candidate_clades[
          i,
          ,
          drop = FALSE
        ]
      )

      selected_nodes <- c(
        selected_nodes,
        current_node
      )

    }

  }

  candidate_clades <- selected_clades

}


# ============================================================
# REPORT GENOTYPE CLADES
# ============================================================

if (nrow(candidate_clades) == 0) {

  cat(
    "No genotype-consistent clades identified.\n"
  )

} else {

  cat(
    "Genotype-consistent clades identified:",
    nrow(candidate_clades),
    "\n\n"
  )

  print(
    candidate_clades
  )

}


# ============================================================
# CREATE GGTREE OBJECT
# ============================================================

p <- ggtree(
  tree
)

p$data$label <- trimws(
  p$data$label
)


# ============================================================
# ADD METADATA TO TREE
# ============================================================

p$data <- left_join(
  p$data,
  metadata,
  by = c(
    "label" = SEQUENCE_ID_COLUMN
  )
)


# ============================================================
# DETERMINE HIGHLIGHT EXTENSION
# ============================================================

tree_x_max <- max(
  p$data$x,
  na.rm = TRUE
)

extend_to <- tree_x_max *
  HIGHLIGHT_EXTEND_FRACTION


# ============================================================
# DEFINE GENOTYPE COLOURS
# ============================================================

if (nrow(candidate_clades) > 0) {

  genotype_names <- unique(
    candidate_clades$genotype
  )

  pastel_base <- c(
    "#A6CEE3",
    "#B2DF8A",
    "#CAB2D6",
    "#FDBF6F",
    "#FFFFB3",
    "#FB9A99",
    "#80CDC1",
    "#D9D9D9",
    "#BC80BD",
    "#CCEBC5"
  )

  if (
    length(genotype_names) <=
    length(pastel_base)
  ) {

    genotype_colours <- pastel_base[
      seq_along(genotype_names)
    ]

  } else {

    genotype_colours <- colorRampPalette(
      pastel_base
    )(
      length(genotype_names)
    )

  }

  names(genotype_colours) <- genotype_names

} else {

  genotype_colours <- character(0)

}


# ============================================================
# ADD GENOTYPE HIGHLIGHTS
# ============================================================

if (nrow(candidate_clades) > 0) {

  for (
    i in seq_len(
      nrow(candidate_clades)
    )
  ) {

    node_i <- candidate_clades$node[i]

    genotype_i <- candidate_clades$genotype[i]

    p <- p +
      geom_hilight(
        node = node_i,
        fill = genotype_colours[
          genotype_i
        ],
        alpha = HIGHLIGHT_ALPHA,
        extendto = extend_to
      )

  }

}


# ============================================================
# REDRAW TREE BRANCHES
# ============================================================

p <- p +
  geom_tree(
    linewidth = TREE_LINEWIDTH,
    colour = TREE_COLOUR
  )


# ============================================================
# TIP COLOUR VARIABLE HANDLING
# ============================================================
#
# IMPORTANT:
# TIP_COLOUR_COLUMN is the ONLY user setting that needs to be changed
# when you want to colour tips by a different metadata column.
#
# Blank/unknown values are converted to "NA" ONLY for visualisation
# and are displayed in grey. The original metadata CSV is NOT changed.
# ============================================================

if (!TIP_COLOUR_COLUMN %in% colnames(p$data)) {

  stop(
    paste(
      "TIP_COLOUR_COLUMN not found in metadata:",
      TIP_COLOUR_COLUMN
    )
  )

}

tip_colour_vector <- as.character(
  p$data[[TIP_COLOUR_COLUMN]]
)

tip_colour_vector <- trimws(
  tip_colour_vector
)

tip_colour_vector[
  is.na(tip_colour_vector) |
    tip_colour_vector == ""
] <- "NA"

p$data[[TIP_COLOUR_COLUMN]] <- tip_colour_vector


# ============================================================
# GET TIP-COLOUR VALUES REPRESENTED AMONG TREE TIPS
# ============================================================

tip_vector <- as.logical(
  p$data$isTip
)

tip_colour_values <- tip_colour_vector[
  tip_vector
]

tip_colour_values <- unique(
  tip_colour_values
)

tip_colour_values <- sort(
  tip_colour_values
)

# Put unknown values at the end.

if ("NA" %in% tip_colour_values) {

  tip_colour_values <- tip_colour_values[
    tip_colour_values != "NA"
  ]

  tip_colour_values <- c(
    tip_colour_values,
    "NA"
  )

}


# ============================================================
# CREATE TIP-COLOUR PALETTE
# ============================================================

if (length(tip_colour_values) > 0) {

  # Use a high-contrast qualitative palette so that a small number
  # of regions remain clearly distinguishable from one another.
  # The palette is intentionally independent of the number of regions
  # beyond selecting the required number of colours.
  distinct_region_palette <- c(
    "#6BAED6",
    "#F4A261",
    "#74C69D",
    "#C8A2C8",
    "#E9C46A",
    "#8EC6E8",
    "#E88B9B",
    "#A8A8A8",
    "#C49A6C",
    "#6A3D9A",
    "#1B9E77",
    "#E31A1C"
  )

  if (length(tip_colour_values) <= length(distinct_region_palette)) {

    tip_colour_palette <- distinct_region_palette[
      seq_along(tip_colour_values)
    ]

  } else {

    tip_colour_palette <- grDevices::colorRampPalette(
      distinct_region_palette
    )(
      length(tip_colour_values)
    )

  }

  names(tip_colour_palette) <- tip_colour_values

  # Unknown value = grey
  if ("NA" %in% names(tip_colour_palette)) {
    tip_colour_palette["NA"] <- "grey85"
  }

} else {

  tip_colour_palette <- c(
    "NA" = "grey85"
  )

}


# ============================================================
# ADD TIP COLOURS
# ============================================================

p <- p +
  geom_tippoint(
    aes(
      colour = .data[[TIP_COLOUR_COLUMN]]
    ),
    size = TIP_SIZE
  ) +
  scale_colour_manual(
    values = tip_colour_palette,
    na.value = "grey85",
    name = TIP_COLOUR_COLUMN,
    drop = FALSE,
    guide = guide_legend(nrow = TIP_LEGEND_NROW)
  )


# ============================================================
# IDENTIFY STUDY SAMPLES
# ============================================================

study_data <- p$data[
  p$data$isTip &
    !is.na(
      p$data[[STUDY_SAMPLE_COLUMN]]
    ) &
    tolower(
      p$data[[STUDY_SAMPLE_COLUMN]]
    ) ==
    tolower(
      STUDY_SAMPLE_VALUE
    ),
  ,
  drop = FALSE
]

cat(
  "Study sample tips identified:",
  nrow(study_data),
  "\n"
)


# ============================================================
# ADD FILLED HIGHLIGHTS TO STUDY SAMPLES
# ============================================================
#
# The normal country-coloured tip is retained underneath.
# Each study sample then receives a slightly larger filled point
# with its own uniquely coloured border.
#
# Shape 21 is used because fill and colour are independently
# controllable.
# ============================================================

if (nrow(study_data) > 0) {

  # Use one consistent bold border colour for ALL study samples.
  # The fill follows the same metadata variable used for tip colouring.
  STUDY_BORDER_COLOUR <- "#D73027"
  study_data$study_border_colour <- STUDY_BORDER_COLOUR

  # Keep the fill consistent with the metadata variable used for tip colours.
  study_tip_colour <- as.character(
    study_data[[TIP_COLOUR_COLUMN]]
  )

  study_data$study_fill_colour <- unname(
    tip_colour_palette[study_tip_colour]
  )

  study_data$study_fill_colour[
    is.na(study_data$study_fill_colour)
  ] <- "grey85"

  p <- p +
    geom_point(
      data = study_data,
      aes(
        x = x,
        y = y
      ),
      shape = STUDY_POINT_SHAPE,
      size = STUDY_POINT_SIZE,
      stroke = STUDY_POINT_STROKE,
      fill = study_data$study_fill_colour,
      colour = study_data$study_border_colour,
      inherit.aes = FALSE
    )

}


# ============================================================
# LABEL GENOTYPE HIGHLIGHTS ON THE TREE
# ============================================================
#
# IMPORTANT:
# The genotype detection logic above is unchanged.
# The labels below are added only after the existing clades have
# been selected, so the three existing genotype highlights remain
# unchanged.
# ============================================================

if (nrow(candidate_clades) > 0) {

  genotype_label_data <- p$data[
    !p$data$isTip &
      p$data$node %in% candidate_clades$node,
    ,
    drop = FALSE
  ]

  genotype_label_data <- left_join(
    genotype_label_data,
    candidate_clades[
      ,
      c(
        "node",
        "genotype"
      ),
      drop = FALSE
    ],
    by = "node"
  )

  genotype_label_data$label_x <- genotype_label_data$x +
    (tree_x_max * 0.012)

  p <- p +
    geom_text(
      data = genotype_label_data,
      aes(
        x = label_x,
        y = y,
        label = genotype
      ),
      hjust = 0,
      vjust = 0.5,
      size = 3.2,
      fontface = "bold",
      inherit.aes = FALSE
    )

}


# ============================================================
# FINAL TREE THEME
# ============================================================

p <- p +
  theme_tree2(
    legend.position = "bottom"
  ) +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5,
      margin = margin(
        b = 18
      )
    ),
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    ),
    legend.key.size = grid::unit(
      0.45,
      "cm"
    ),
    plot.margin = margin(
      10,
      10,
      10,
      10
    )
  ) +
  guides(
    colour = guide_legend(
      nrow = TIP_LEGEND_NROW,
      byrow = TRUE
    )
  ) +
  ggtitle(
    "Maximum-likelihood phylogeny of PPR virus whole-genome sequences"
  )


# ============================================================
# SAVE FIGURE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SAVING FIGURE\n")
cat("============================================================\n")

ggsave(
  filename = OUTPUT_FILE,
  plot = p,
  width = OUTPUT_WIDTH,
  height = OUTPUT_HEIGHT,
  units = "in",
  dpi = OUTPUT_DPI,
  device = "png"
)


# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat(
  "Analysis:",
  analysis_name,
  "\n"
)

cat(
  "Genotype mode:",
  GENOTYPE_MODE,
  "\n"
)

cat(
  "Final tree tips:",
  length(tree$tip.label),
  "\n"
)

cat(
  "Genotype-highlighted clades:",
  nrow(candidate_clades),
  "\n"
)

cat(
  "Study samples:",
  nrow(study_data),
  "\n"
)

cat(
  "Output file:",
  OUTPUT_FILE,
  "\n"
)

cat("============================================================\n")