# Load count data
counts <- read.csv("GSE248260_white_matter_count_matrix.csv", 
                   header = TRUE, row.names = 1, check.names = FALSE)
# Load biomaRt
#if (!requireNamespace("biomaRt", quietly = TRUE)) {
  #install.packages("BiocManager")
  #BiocManager::install("biomaRt")
#}
library(biomaRt)

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Clean rownames (remove version numbers like .10 from ENSG000001234.10)
gene_ids <- sub("\\..*", "", rownames(counts))

# Get gene lengths
gene_info <- getBM(attributes = c("ensembl_gene_id", "transcript_length"),
                   filters = "ensembl_gene_id",
                   values = gene_ids,
                   mart = ensembl)

# Use median transcript length per gene
library(dplyr)
gene_info <- getBM(attributes = c("ensembl_gene_id", "transcript_length"),
                   filters = "ensembl_gene_id",
                   values = gene_ids,
                   mart = ensembl)
# Calculate median transcript length for each gene
library(dplyr)

gene_lengths <- gene_info %>%
  group_by(ensembl_gene_id) %>%
  summarise(length = median(transcript_length, na.rm = TRUE))

# Match and align with counts matrix
counts$ensembl_gene_id <- sub("\\..*", "", rownames(counts))
tpm_input <- merge(counts, gene_lengths, by = "ensembl_gene_id")
rownames(tpm_input) <- tpm_input$ensembl_gene_id
tpm_input <- tpm_input[ , !(colnames(tpm_input) %in% c("ensembl_gene_id"))]
calculateTPM <- function(counts, lengths) {
  rpk <- counts / (lengths / 1000)
  scaling_factors <- colSums(rpk, na.rm = TRUE) / 1e6
  tpm <- sweep(rpk, 2, scaling_factors, "/")
  return(tpm)
}

# Run TPM normalization
lengths <- tpm_input$length
counts_only <- tpm_input[, -which(colnames(tpm_input) == "length")]
tpm_matrix <- calculateTPM(counts_only, lengths)

# Save TPM output
write.csv(tpm_matrix, "GSE248260_TPM_normalized.csv")
