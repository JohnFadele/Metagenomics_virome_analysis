# ============================================================
# STEP 5: EXTRACT ALIGNMENT-SPECIFIC METADATA
#
# Match sequence IDs from the trimmed FASTA against the
# sequence_id column in the general metadata file.
#
# ============================================================


# ============================================================
# USER SETTINGS & INPUT DEFINITIONS
# ============================================================

# Base project directory
PROJECT_DIR <- "/home/fadele/ARS_NACA/ppr_phylo"

# 1. INPUT METADATA FILE (Specify explicit path)
INPUT_METADATA_FILE <- file.path(
  PROJECT_DIR,
  "metadata.csv"
)

# 2. INPUT ALIGNMENT FASTA (Specify explicit path)
INPUT_ALIGNMENT_FILE <- file.path(
  PROJECT_DIR,
  "02_alignment",
  "ppr_pub_study_og_trimmed.fasta"
)

# 3. METADATA COLUMN SETTINGS
METADATA_ID_COLUMN <- "sequence_id"

# 4. OUTPUT DIRECTORY
OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "04_metadata"
)


# ============================================================
# CREATE OUTPUT DIRECTORY & LOG FILE SETUP
# ============================================================

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

LOG_FILE <- file.path(
  OUTPUT_DIR,
  "metadata_extraction.log"
)

if (file.exists(LOG_FILE)) {
  file.remove(LOG_FILE)
}

log_message <- function(...) {
  message_text <- paste0(...)
  cat(message_text, "\n", file = LOG_FILE, append = TRUE)
  cat(message_text, "\n")
}


# ============================================================
# INPUT VALIDATION
# ============================================================

# Verify metadata file exists
if (!file.exists(INPUT_METADATA_FILE)) {
  stop(
    paste0(
      "\n[ERROR] Metadata file not found at specified location:\n",
      INPUT_METADATA_FILE,
      "\n"
    )
  )
}

# Verify alignment FASTA exists
if (!file.exists(INPUT_ALIGNMENT_FILE)) {
  stop(
    paste0(
      "\n[ERROR] Alignment FASTA file not found at specified location:\n",
      INPUT_ALIGNMENT_FILE,
      "\n"
    )
  )
}


# ============================================================
# DETERMINE ANALYSIS NAME & DEFINE OUTPUT FILES
# ============================================================

# Derive analysis prefix from the FASTA file name
fasta_filename <- basename(INPUT_ALIGNMENT_FILE)
analysis_name <- sub("_trimmed\\.fasta$", "", fasta_filename)

# Fallback if filename does not contain '_trimmed.fasta'
if (analysis_name == fasta_filename) {
  analysis_name <- tools::file_path_sans_ext(fasta_filename)
}

FILTERED_METADATA_FILE <- file.path(
  OUTPUT_DIR,
  paste0(analysis_name, "_alignment_metadata.csv")
)

MISSING_METADATA_FILE <- file.path(
  OUTPUT_DIR,
  paste0(analysis_name, "_metadata_missing.csv")
)

UNUSED_METADATA_FILE <- file.path(
  OUTPUT_DIR,
  paste0(analysis_name, "_metadata_unused.csv")
)


# ============================================================
# START LOG
# ============================================================

log_message("============================================================")
log_message("ALIGNMENT-SPECIFIC METADATA EXTRACTION")
log_message("============================================================")
log_message("Start time: ", Sys.time())
log_message("Project directory: ", PROJECT_DIR)
log_message("Input metadata file: ", INPUT_METADATA_FILE)
log_message("Input alignment FASTA: ", INPUT_ALIGNMENT_FILE)
log_message("Metadata ID column: ", METADATA_ID_COLUMN)
log_message("Output directory: ", OUTPUT_DIR)
log_message("")


# ============================================================
# READ ALIGNMENT FASTA & EXTRACT IDs
# ============================================================

log_message("Reading alignment FASTA...")

fasta_lines <- readLines(
  INPUT_ALIGNMENT_FILE,
  warn = FALSE
)

# Extract FASTA headers
headers <- fasta_lines[grepl("^>", fasta_lines)]

if (length(headers) == 0) {
  stop("ERROR: No FASTA headers were found in the alignment file.")
}

# Clean and extract sequence IDs from headers
seq_ids <- sub("^>", "", headers)
seq_ids <- gsub("\r", "", seq_ids)
seq_ids <- sapply(strsplit(seq_ids, "\\|"), `[`, 1)
seq_ids <- gsub("^_R_", "", seq_ids)
seq_ids <- trimws(seq_ids)

# Check for empty sequence IDs
empty_ids <- which(is.na(seq_ids) | seq_ids == "")
if (length(empty_ids) > 0) {
  stop(
    paste0(
      "ERROR: ",
      length(empty_ids),
      " FASTA header(s) produced an empty sequence ID."
    )
  )
}

# Check for duplicate IDs in alignment
duplicate_ids <- unique(seq_ids[duplicated(seq_ids)])
if (length(duplicate_ids) > 0) {
  log_message("WARNING: Duplicate sequence IDs detected in alignment:")
  log_message(paste(duplicate_ids, collapse = ", "))
  log_message("")
}


# ============================================================
# READ & VALIDATE METADATA
# ============================================================

log_message("Reading metadata...")

metadata <- read.csv(
  INPUT_METADATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!METADATA_ID_COLUMN %in% names(metadata)) {
  stop(
    paste0(
      "ERROR: The metadata file does not contain the column '",
      METADATA_ID_COLUMN,
      "'.\n\n",
      "Available columns are:\n",
      paste(names(metadata), collapse = "\n")
    )
  )
}

# Clean sequence IDs in metadata
metadata[[METADATA_ID_COLUMN]] <- trimws(
  as.character(metadata[[METADATA_ID_COLUMN]])
)

# Check for duplicate metadata sequence_id values
duplicate_metadata_ids <- unique(
  metadata[[METADATA_ID_COLUMN]][duplicated(metadata[[METADATA_ID_COLUMN]])]
)

if (length(duplicate_metadata_ids) > 0) {
  log_message("WARNING: Duplicate sequence_id values detected in metadata:")
  log_message(paste(duplicate_metadata_ids, collapse = ", "))
  log_message("The first matching metadata record will be retained.")
  log_message("")
}


# ============================================================
# MATCH FASTA IDs TO METADATA & GENERATE TABLES
# ============================================================

metadata_index <- match(
  seq_ids,
  metadata[[METADATA_ID_COLUMN]]
)

# Identify missing metadata
missing_indices <- which(is.na(metadata_index))
missing_ids <- seq_ids[missing_indices]

# Retain matched metadata preserving FASTA sequence order
matched_indices <- which(!is.na(metadata_index))
filtered_metadata <- metadata[metadata_index[matched_indices], , drop = FALSE]

matched_seq_ids <- seq_ids[matched_indices]
filtered_metadata[[METADATA_ID_COLUMN]] <- matched_seq_ids
rownames(filtered_metadata) <- NULL

# Create missing-metadata data frame
if (length(missing_ids) > 0) {
  missing_metadata <- data.frame(
    sequence_id = missing_ids,
    stringsAsFactors = FALSE
  )
} else {
  missing_metadata <- data.frame(
    sequence_id = character(0),
    stringsAsFactors = FALSE
  )
}

# Identify unused metadata
used_ids <- unique(filtered_metadata[[METADATA_ID_COLUMN]])
unused_metadata <- metadata[
  !(metadata[[METADATA_ID_COLUMN]] %in% used_ids),
  ,
  drop = FALSE
]


# ============================================================
# WRITE OUTPUT FILES
# ============================================================

write.csv(
  filtered_metadata,
  FILTERED_METADATA_FILE,
  row.names = FALSE
)

write.csv(
  missing_metadata,
  MISSING_METADATA_FILE,
  row.names = FALSE
)

write.csv(
  unused_metadata,
  UNUSED_METADATA_FILE,
  row.names = FALSE
)


# ============================================================
# LOG SUMMARY STATISTICS
# ============================================================

total_alignment_sequences <- length(seq_ids)
total_metadata_records <- nrow(metadata)
matched_sequences <- nrow(filtered_metadata)
missing_sequence_count <- length(missing_ids)
unused_metadata_count <- nrow(unused_metadata)

log_message("")
log_message("============================================================")
log_message("SUMMARY")
log_message("============================================================")
log_message("Total sequences in alignment: ", total_alignment_sequences)
log_message("Total metadata records:       ", total_metadata_records)
log_message("Metadata records retained:    ", matched_sequences)
log_message("Alignment sequences missing:  ", missing_sequence_count)
log_message("Metadata records not used:    ", unused_metadata_count)
log_message("")

if (missing_sequence_count > 0) {
  log_message("WARNING: The following alignment sequences have no metadata:")
  for (id in missing_ids) {
    log_message("  ", id)
  }
  log_message("")
} else {
  log_message("All alignment sequences have corresponding metadata.")
  log_message("")
}

log_message("Alignment-specific metadata: ", FILTERED_METADATA_FILE)
log_message("Sequences missing metadata:  ", MISSING_METADATA_FILE)
log_message("Metadata records not used:   ", UNUSED_METADATA_FILE)
log_message("Log file:                    ", LOG_FILE)
log_message("")
log_message("End time: ", Sys.time())
log_message("============================================================")


# ============================================================
# CONSOLE SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("METADATA EXTRACTION COMPLETED\n")
cat("============================================================\n")
cat("Alignment sequences:     ", total_alignment_sequences, "\n")
cat("Metadata records:        ", total_metadata_records, "\n")
cat("Records retained:        ", matched_sequences, "\n")
cat("Missing metadata:        ", missing_sequence_count, "\n")
cat("Unused metadata records: ", unused_metadata_count, "\n\n")

cat("Main output:\n", FILTERED_METADATA_FILE, "\n\n")
cat("Missing metadata report:\n", MISSING_METADATA_FILE, "\n\n")
cat("Unused metadata report:\n", UNUSED_METADATA_FILE, "\n\n")
cat("Log file:\n", LOG_FILE, "\n")
cat("============================================================\n")