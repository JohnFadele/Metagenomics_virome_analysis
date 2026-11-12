# ============================================================
# FIND THE CLOSEST SEQUENCES TO EACH OF FOCAL SEQUENCES
# USING PHYLOGENETIC DISTANCE
# ============================================================


# ============================================================
# USER SETTINGS
# ============================================================

PROJECT_DIR <- "/home/fadele/ARS_NACA/PPR_phylo"

TREE_FILE <- file.path(
  PROJECT_DIR,
  "03_IQTREE",
  "PPR_WG_phylogeny_seq_rooted.treefile"
)

OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "05_visualisation"
)

# Number of closest sequences wanted for each focal sequence
N_CLOSEST <- 20


# ============================================================
# FOCAL SEQUENCES
# ============================================================

FOCAL_SEQUENCES <- c(
  "AD_YS_KR_GT_006",
  "AD_YS_KR_GT_007",
  "AD_YS_KR_SP_013",
  "KN_DL_AB_GT_007_S49",
  "KN_DL_KM_AB_SP_004_S57",
  "KN_DL_KM_AB_SP_007_S59",
  "SK_SN_GK_KR_GT_032",
  "SK_SN_GK_KR_SP_017"
)


# ============================================================
# OUTPUT FILES
# ============================================================

OUTPUT_ALL <- file.path(
  OUTPUT_DIR,
  "20_closest_sequences_per_focal.csv"
)

OUTPUT_UNIQUE <- file.path(
  OUTPUT_DIR,
  "unique_closest_sequences_across_20_focal.csv"
)


# ============================================================
# LOAD PACKAGE
# ============================================================

library(ape)


# ============================================================
# CHECK FILES
# ============================================================

if (!file.exists(TREE_FILE)) {

  stop(
    paste(
      "Tree file not found:",
      TREE_FILE
    )
  )

}

if (!dir.exists(OUTPUT_DIR)) {

  dir.create(
    OUTPUT_DIR,
    recursive = TRUE
  )

}


# ============================================================
# READ TREE
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
  "Number of tree tips:",
  length(tree$tip.label),
  "\n"
)


# ============================================================
# CHECK FOCAL SEQUENCES
# ============================================================

focal_present <- FOCAL_SEQUENCES[
  FOCAL_SEQUENCES %in% tree$tip.label
]

focal_missing <- FOCAL_SEQUENCES[
  !FOCAL_SEQUENCES %in% tree$tip.label
]


if (length(focal_missing) > 0) {

  cat("\n")
  cat("WARNING: The following focal sequences were NOT found:\n")

  print(focal_missing)

}


if (length(focal_present) == 0) {

  stop(
    "None of the focal sequences were found in the tree."
  )

}


cat("\n")
cat(
  "Focal sequences found:",
  length(focal_present),
  "of",
  length(FOCAL_SEQUENCES),
  "\n"
)


# ============================================================
# CALCULATE PAIRWISE PHYLOGENETIC DISTANCES
# ============================================================

cat("\n")
cat("============================================================\n")
cat("CALCULATING PHYLOGENETIC DISTANCES\n")
cat("============================================================\n")

distance_matrix <- cophenetic.phylo(tree)


# ============================================================
# FIND N CLOSEST SEQUENCES FOR EACH FOCAL SEQUENCE
# ============================================================

results <- data.frame(
  focal_sequence = character(),
  rank = integer(),
  closest_sequence = character(),
  phylogenetic_distance = numeric(),
  stringsAsFactors = FALSE
)


for (focal in focal_present) {

  cat(
    "\nFinding closest sequences to:",
    focal,
    "\n"
  )

  # Extract distances from focal sequence
  distances <- distance_matrix[
    focal,
    ]

  # Remove focal sequence itself
  distances <- distances[
    names(distances) != focal
  ]

  # Sort from closest to most distant
  distances <- sort(
    distances,
    decreasing = FALSE
  )

  # Select the requested number
  n_to_take <- min(
    N_CLOSEST,
    length(distances)
  )

  closest <- distances[
    seq_len(n_to_take)
  ]

  focal_results <- data.frame(
    focal_sequence = focal,
    rank = seq_along(closest),
    closest_sequence = names(closest),
    phylogenetic_distance = as.numeric(closest),
    stringsAsFactors = FALSE
  )

  results <- rbind(
    results,
    focal_results
  )

}


# ============================================================
# SAVE COMPLETE RESULTS
# ============================================================

write.csv(
  results,
  OUTPUT_ALL,
  row.names = FALSE
)


# ============================================================
# CREATE UNIQUE LIST OF CLOSEST SEQUENCES
# ============================================================

unique_sequences <- unique(
  results$closest_sequence
)

unique_results <- data.frame(
  sequence_id = sort(unique_sequences),
  stringsAsFactors = FALSE
)

write.csv(
  unique_results,
  OUTPUT_UNIQUE,
  row.names = FALSE
)


# ============================================================
# PRINT RESULTS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("RESULTS\n")
cat("============================================================\n")

print(
  results
)


cat("\n")
cat(
  "Unique sequences identified across all focal sequences:",
  nrow(unique_results),
  "\n"
)

cat("\n")
cat("Output files:\n")

cat(
  "1. ",
  OUTPUT_ALL,
  "\n",
  sep = ""
)

cat(
  "2. ",
  OUTPUT_UNIQUE,
  "\n",
  sep = ""
)


cat("\n")
cat("============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")
