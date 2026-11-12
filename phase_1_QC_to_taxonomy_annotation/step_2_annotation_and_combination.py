#!/usr/bin/env python3
"""
==============================================================================
PHASE 2: BRACKEN TAXONOMIC ANNOTATION & COMBINATION PIPELINE
==============================================================================
Steps:
  1. Annotate Bracken species outputs using NCBI taxonomy (nodes.dmp / names.dmp).
  2. Filter host contaminants (e.g., Homo sapiens) per sample.
  3. Recalculate relative abundance fractions based on non-host reads.
  4. Merge all annotated sample outputs into a single consolidated CSV dataset[cite: 6].
==============================================================================
"""

import os
import csv
from functools import lru_cache

# ==============================================================================
# 1. USER CONFIGURATION
# Set your working directories, taxonomy files, and filtering settings here.
# ==============================================================================

# Database containing NCBI nodes.dmp and names.dmp
TAXONOMY_DIR = "/home/fadele/Kraken_full/k2_standard_20260626"

# Directory containing step 4 Bracken output (*.species.bracken.txt)
BRACKEN_DIR = "/home/fadele/ppr_extra/results/bracken"

# Directory where per-sample annotated TSVs will be saved
ANNOTATED_OUTPUT_DIR = "/home/fadele/ppr_extra/results/annotated_tsv"

# Path for the final combined abundance matrix
FINAL_COMBINED_CSV = "/home/fadele/ppr_extra/results/combined_species_abundance.csv"

# Host Filter Settings
FILTER_HUMAN = True
HUMAN_TAXON_NAME = "Homo sapiens"

# Major Taxonomic Domain Classification
TOP_GROUPS = {"Viruses", "Bacteria", "Archaea", "Eukaryota"}

# ==============================================================================
# 2. SETUP & DIRECTORY CREATION
# ==============================================================================

os.makedirs(ANNOTATED_OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.dirname(FINAL_COMBINED_CSV), exist_ok=True)

NODES_FILE = os.path.join(TAXONOMY_DIR, "nodes.dmp")
NAMES_FILE = os.path.join(TAXONOMY_DIR, "names.dmp")

# Check required database files exist
if not os.path.isfile(NODES_FILE) or not os.path.isfile(NAMES_FILE):
    raise FileNotFoundError(
        f"[ERROR] Could not find 'nodes.dmp' or 'names.dmp' in: {TAXONOMY_DIR}"
    )

print("==============================================================================")
print("STARTING PHASE 2: TAXONOMIC ANNOTATION & COMBINATION")
print(f"Bracken Directory : {BRACKEN_DIR}")
print(f"Annotated Directory: {ANNOTATED_OUTPUT_DIR}")
print(f"Final CSV Output   : {FINAL_COMBINED_CSV}")
print("==============================================================================")

# ==============================================================================
# 3. LOAD NCBI TAXONOMY
# ==============================================================================

print("\nLoading NCBI taxonomy dump files...")

parent_map = {}
rank_map = {}
name_map = {}

with open(NODES_FILE, "r") as f:
    for line in f:
        fields = [x.strip() for x in line.split("|")]
        taxid = int(fields[0])
        parent_map[taxid] = int(fields[1])
        rank_map[taxid] = fields[2]

with open(NAMES_FILE, "r") as f:
    for line in f:
        fields = [x.strip() for x in line.split("|")]
        if fields[3] == "scientific name":
            name_map[int(fields[0])] = fields[1]

print(f"[DONE] Loaded {len(name_map):,} scientific names.")

# ==============================================================================
# 4. TAXONOMIC LINEAGE TRAVERSAL
# ==============================================================================

@lru_cache(maxsize=None)
def get_annotation(taxid):
    """Traverses the taxonomy tree upwards to locate Family, Species, and Top Domain."""
    family = ""
    species = ""
    highest = ""

    current = taxid
    while True:
        current_name = name_map.get(current, "")
        current_rank = rank_map.get(current, "")

        if current_rank == "family" and not family:
            family = current_name

        if current_rank == "species" and not species:
            species = current_name

        if current_name in TOP_GROUPS:
            highest = current_name

        if current == 1 or current not in parent_map:
            break

        current = parent_map.get(current, 1)

    return highest, family, species

# ==============================================================================
# 5. STEP 5 & 6 UNIFIED PROCESSING
# ==============================================================================

bracken_files = sorted([
    f for f in os.listdir(BRACKEN_DIR) 
    if f.endswith(".species.bracken.txt")
])

if not bracken_files:
    raise FileNotFoundError(f"No '.species.bracken.txt' files found in: {BRACKEN_DIR}")

print(f"\nFound {len(bracken_files)} Bracken reports to annotate and combine.")

all_combined_rows = []
total_combined_rows = 0
removed_human_rows = 0
removed_human_reads = 0.0

fieldnames_tsv = [
    "name",
    "taxonomy_id",
    "taxonomy_lvl",
    "kraken_assigned_reads",
    "added_reads",
    "new_est_reads",
    "fraction_total_reads",
    "Highest_taxon",
    "Family",
    "Species"
]

fieldnames_csv = ["Sample_ID"] + fieldnames_tsv

for file in bracken_files:
    # Deriving Sample Name
    sample_id = file.replace(".species.bracken.txt", "")
    sample_id = sample_id.replace("_clean", "")

    infile = os.path.join(BRACKEN_DIR, file)
    annotated_tsv = os.path.join(ANNOTATED_OUTPUT_DIR, f"{sample_id}.annotated.tsv")

    print(f"[PROCESSING] {sample_id}")

    sample_rows = []

    # --------------------------------------------------------------------------
    # Step A: Read Bracken File & Annotate
    # --------------------------------------------------------------------------
    with open(infile, "r") as fin:
        reader = csv.DictReader(fin, delimiter="\t")

        for row in reader:
            try:
                taxid = int(row["taxonomy_id"])
            except (ValueError, KeyError):
                continue

            highest, family, species = get_annotation(taxid)

            row_dict = {
                "name": row.get("name", ""),
                "taxonomy_id": row.get("taxonomy_id", ""),
                "taxonomy_lvl": row.get("taxonomy_lvl", ""),
                "kraken_assigned_reads": row.get("kraken_assigned_reads", ""),
                "added_reads": row.get("added_reads", ""),
                "new_est_reads": row.get("new_est_reads", ""),
                "fraction_total_reads": row.get("fraction_total_reads", ""),
                "Highest_taxon": highest,
                "Family": family,
                "Species": species
            }

            # ------------------------------------------------------------------
            # Step B: Filter Host / Human Reads
            # ------------------------------------------------------------------
            if FILTER_HUMAN:
                is_human = (
                    row_dict["Species"].strip() == HUMAN_TAXON_NAME or
                    row_dict["name"].strip() == HUMAN_TAXON_NAME
                )
                if is_human:
                    removed_human_rows += 1
                    try:
                        removed_human_reads += float(row_dict["new_est_reads"])
                    except (ValueError, TypeError):
                        pass
                    continue

            sample_rows.append(row_dict)

    # --------------------------------------------------------------------------
    # Step C: Write Per-Sample Annotated TSV File
    # --------------------------------------------------------------------------
    with open(annotated_tsv, "w", newline="") as fout_tsv:
        writer_tsv = csv.DictWriter(fout_tsv, fieldnames=fieldnames_tsv, delimiter="\t")
        writer_tsv.writeheader()
        writer_tsv.writerows(sample_rows)

    # --------------------------------------------------------------------------
    # Step D: Recalculate Relative Abundance Fractions (Non-Host Reads)
    # --------------------------------------------------------------------------
    nonhost_total_reads = sum(
        float(r["new_est_reads"]) for r in sample_rows 
        if r["new_est_reads"] not in ("", None)
    )

    for r in sample_rows:
        try:
            est_reads = float(r["new_est_reads"])
        except (ValueError, TypeError):
            est_reads = 0.0

        if nonhost_total_reads > 0:
            new_fraction = est_reads / nonhost_total_reads
        else:
            new_fraction = 0.0

        r["fraction_total_reads"] = f"{new_fraction:.8f}"
        
        # Prepare row for combined CSV
        combined_row = {"Sample_ID": sample_id}
        combined_row.update(r)
        all_combined_rows.append(combined_row)
        total_combined_rows += 1

# ==============================================================================
# 6. WRITE FINAL COMBINED MATRIX CSV
# ==============================================================================

with open(FINAL_COMBINED_CSV, "w", newline="") as fout_csv:
    writer_csv = csv.DictWriter(fout_csv, fieldnames=fieldnames_csv)
    writer_csv.writeheader()
    writer_csv.writerows(all_combined_rows)

print("\n==============================================================================")
print("PHASE 2 PIPELINE COMPLETE")
print(f"Total Samples Processed : {len(bracken_files)}")
print(f"Total Species Entries   : {total_combined_rows}")
if FILTER_HUMAN:
    print(f"Human Taxa Removed      : {removed_human_rows} entries ({removed_human_reads:,.0f} reads)")
print(f"Per-Sample TSVs Saved   : {ANNOTATED_OUTPUT_DIR}")
print(f"Combined Abundance CSV  : {FINAL_COMBINED_CSV}")
print("==============================================================================")