# ===================================================================
# QC Script for GenomeStudio SNP Array Data
# Author: James Standish
# Purpose: Compute marker- and individual-level QC metrics for
#          genotype data using phenotype mapping information.
# ===================================================================

# -----------------------------
# Load packages
# -----------------------------
library(tidyverse)
library(ggplot2)

# -----------------------------
# 1. File paths
# -----------------------------
geno_file  <- "C:/Users/StandishJ/Downloads/FinalProject/Rqtl_Inputs/R5xRG_Genome_Studio_Carolina_08_2016.txt"
pheno_file <- "C:/Users/StandishJ/Downloads/FinalProject/Rqtl_Inputs/Phenotype_Genotype_Mapping.txt"
output_dir <- "C:/Users/StandishJ/Downloads/FinalProject/QC_Outputs"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. Read files
# -----------------------------
geno  <- read.delim(geno_file, stringsAsFactors = FALSE, check.names = FALSE)
pheno <- read.delim(pheno_file, stringsAsFactors = FALSE, check.names = FALSE)

# -----------------------------
# 3. Validate expected columns
# -----------------------------
required_pheno_cols <- c("Genotype", "MeanDiseaseIncidence", "TotalFlowerClusters", "CleanID")
missing_pheno_cols <- setdiff(required_pheno_cols, names(pheno))
if (length(missing_pheno_cols) > 0) {
  stop("Phenotype file is missing required columns: ", paste(missing_pheno_cols, collapse = ", "))
}

required_geno_cols <- c("Index", "Name", "Chr", "Position")
missing_geno_cols <- setdiff(required_geno_cols, names(geno))
if (length(missing_geno_cols) > 0) {
  stop("Genotype file is missing required columns: ", paste(missing_geno_cols, collapse = ", "))
}

# -----------------------------
# 4. Clean genotype column names
# -----------------------------
if (ncol(geno) < 10) {
  stop("Genotype file does not appear to have the expected structure.")
}

geno_cols <- colnames(geno)[5:(ncol(geno) - 5)]  # genotype/sample columns
geno_cols_clean <- geno_cols %>%
  trimws() %>%
  gsub("\\s*\\[.*\\]", "", .) %>%  # remove [number]
  gsub("\\.", "", .)               # remove dots

colnames(geno)[5:(ncol(geno) - 5)] <- geno_cols_clean

# -----------------------------
# 5. Clean phenotype IDs
# -----------------------------
pheno <- pheno %>%
  mutate(
    CleanID = as.character(CleanID) %>%
      trimws() %>%
      gsub("\\s*\\[.*\\]", "", .) %>%
      gsub("\\.", "", .),
    Genotype = as.character(Genotype)
  )

# -----------------------------
# 6. QC: check matching IDs
# -----------------------------
matching_ids <- intersect(pheno$CleanID, geno_cols_clean)
missing_ids   <- setdiff(pheno$CleanID, geno_cols_clean)

if (length(missing_ids) > 0) {
  warning("These phenotype IDs were missing in genotype file: ", paste(missing_ids, collapse = ", "))
}

if (length(matching_ids) == 0) {
  stop("No matching genotype IDs were found between phenotype and genotype files.")
}

# -----------------------------
# 7. Filter genotype to matching IDs
# -----------------------------
geno_filtered <- geno %>%
  select(Index, Name, Chr, Position, all_of(matching_ids))

# -----------------------------
# 8. Convert "NC" to NA
# -----------------------------
geno_filtered <- geno_filtered %>%
  mutate(across(all_of(matching_ids), ~ na_if(trimws(as.character(.)), "NC")))

# -----------------------------
# 9. QC: calculate call rates
# -----------------------------
# Sample call rate
sample_callrate <- colSums(!is.na(geno_filtered[, matching_ids, drop = FALSE])) / nrow(geno_filtered)
sample_callrate_df <- data.frame(
  Sample = names(sample_callrate),
  CallRate = as.numeric(sample_callrate)
)

write.table(
  sample_callrate_df,
  file = file.path(output_dir, "Sample_CallRates.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# SNP call rate
snp_callrate <- rowSums(!is.na(geno_filtered[, matching_ids, drop = FALSE])) / length(matching_ids)
snp_callrate_df <- geno_filtered %>%
  select(Index, Name, Chr, Position) %>%
  mutate(CallRate = as.numeric(snp_callrate))

write.table(
  snp_callrate_df,
  file = file.path(output_dir, "SNP_CallRates.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# -----------------------------
# 10. Plot QC distributions
# -----------------------------
plot_histogram <- function(df, col, title, filename) {
  p <- ggplot(df, aes(x = .data[[col]])) +
    geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
    theme_minimal() +
    labs(title = title, x = "Call Rate", y = "Count")
  
  ggsave(file.path(output_dir, paste0(filename, ".png")), plot = p, width = 6, height = 4)
  ggsave(file.path(output_dir, paste0(filename, ".pdf")), plot = p, width = 6, height = 4)
}

plot_histogram(sample_callrate_df, "CallRate", "Sample Call Rate Distribution", "Sample_CallRate_Histogram")
plot_histogram(snp_callrate_df, "CallRate", "SNP Call Rate Distribution", "SNP_CallRate_Histogram")

# -----------------------------
# 11. Merge phenotype info
# -----------------------------
geno_summary <- pheno %>%
  filter(CleanID %in% matching_ids) %>%
  arrange(CleanID) %>%
  select(Genotype, MeanDiseaseIncidence, TotalFlowerClusters, CleanID)

write.table(
  geno_summary,
  file = file.path(output_dir, "Phenotype_Genotype_Mapping_Cleaned.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# -----------------------------
# 12. Summary report
# -----------------------------
qc_summary <- list(
  Total_Samples = length(matching_ids),
  Missing_Samples = length(missing_ids),
  Total_SNPs = nrow(geno_filtered),
  Sample_CallRate_Mean = mean(sample_callrate, na.rm = TRUE),
  Sample_CallRate_Min = min(sample_callrate, na.rm = TRUE),
  Sample_CallRate_Max = max(sample_callrate, na.rm = TRUE),
  SNP_CallRate_Mean = mean(snp_callrate, na.rm = TRUE),
  SNP_CallRate_Min = min(snp_callrate, na.rm = TRUE),
  SNP_CallRate_Max = max(snp_callrate, na.rm = TRUE)
)

qc_summary_df <- as.data.frame(t(as.data.frame(qc_summary)))

write.table(
  qc_summary_df,
  file = file.path(output_dir, "QC_Summary.txt"),
  sep = "\t", row.names = TRUE, quote = FALSE
)

# -----------------------------
# 13. Genotype distribution summary
# -----------------------------
geno_matrix <- geno_filtered[, matching_ids, drop = FALSE]

geno_counts <- sapply(c("AA", "AB", "BB", "NC"), function(gt) {
  if (gt == "NC") {
    sum(is.na(geno_matrix))
  } else {
    sum(geno_matrix == gt, na.rm = TRUE)
  }
})

total_calls <- nrow(geno_matrix) * ncol(geno_matrix)
geno_prop <- geno_counts / total_calls * 100

geno_summary_table <- data.frame(
  Genotype = names(geno_counts),
  Count = as.numeric(geno_counts),
  Proportion = round(as.numeric(geno_prop), 2)
)

write.table(
  geno_summary_table,
  file = file.path(output_dir, "Genotype_Distribution.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

print(geno_summary_table)

cat("QC pipeline completed. Outputs saved to: ", output_dir, "\n")
cat("Phenotype file CleanID column used successfully.\n")