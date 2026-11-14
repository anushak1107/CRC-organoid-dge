# ==============================================================================
# Script 01: Data Preprocessing & Automated Metadata Parsing (GSE349730)
# ==============================================================================

# Create output directories if they do not exist
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

# 1. Load raw count matrix
raw_counts <- read.csv("data/349730dge.csv", row.names = 1, check.names = FALSE)
cat("Loaded count matrix:", nrow(raw_counts), "genes x", ncol(raw_counts), "samples\n")

# 2. Programmatically parse sample metadata from column names
# Example column format: "P 23 CRC_MG_1" or "P 23 CRC_MGCOL_5FU_1"
sample_names <- colnames(raw_counts)

coldata <- data.frame(
  row.names = sample_names,
  Patient   = sub("^(P \\d+).*", "\\1", sample_names),
  Matrix    = ifelse(grepl("MGCOL", sample_names), "MGCOL", "MG"),
  Treatment = ifelse(grepl("5FU", sample_names), "5FU", "Control"),
  stringsAsFactors = TRUE
)

# Set MG (Matrigel) as the reference baseline level
coldata$Matrix <- relevel(coldata$Matrix, ref = "MG")

# 3. Filter low-count genes (retain genes with >= 10 counts in >= 3 samples)
keep <- rowSums(raw_counts >= 10) >= 3
filtered_counts <- raw_counts[keep, ]
cat("Genes retained after filtering:", nrow(filtered_counts), "\n")

# Save preprocessed objects for downstream scripts
saveRDS(filtered_counts, "results/filtered_counts.rds")
saveRDS(coldata, "results/coldata.rds")
cat("Preprocessing complete! Saved RDS objects to results/\n")
