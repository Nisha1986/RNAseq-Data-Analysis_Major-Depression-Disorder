library(dplyr)
library(tidyverse)
library(GEOquery)
library(GenomicFeatures)
library(txdbmaker)
library(AnnotationDbi)
library(rtracklayer)
#read the data
dat <- read.csv("GSE248260_TPM_normalized.csv",header = TRUE, row.names = 1, check.names = FALSE)
dim(dat)
#get metadata
gse<-getGEO(GEO="GSE248260",GSEMatrix = TRUE)
Sys.setenv("VROOM_CONNECTION_SIZE"= 131072 * 1000)
gse
metadata <-pData(phenoData(gse[[1]]))
head(metadata)
write.csv(metadata, "metadata.csv", row.names = FALSE)
metadata.subset <- metadata[, c(1, 10, 11, 45)]
# Check column names
colnames(metadata.subset)
# Rename columns
colnames(metadata.subset)[colnames(metadata.subset) == "title"] <- "sample"  # GSM ID to Sample
colnames(metadata.subset)[colnames(metadata.subset) == "characteristics_ch1"] <- "tissue"
colnames(metadata.subset)[colnames(metadata.subset) == "characteristics_ch1.1"] <- "diagnosis"
colnames(metadata.subset)[colnames(metadata.subset) == "suicide:ch1"] <- "suicide"
# Print updated column names to verify
print(colnames(metadata.subset))
# Modify the dataset using mutate
metadata.modified <- metadata.subset %>%
  mutate(
    tissue = gsub("^tissue: ", "", tissue),  # Removes "tissue: " from tissue column
    diagnosis = gsub("^diagnosis: ", "", diagnosis),  # Removes "diagnosis: " from diagnosis column
    suicide = gsub("^suicide: ", "", suicide)  # Removes "suicide: " from suicide column
  )
# View the modified dataset
head(metadata.modified)
# Ensure data frame
dat <- as.data.frame(dat)

# Rename GeneName column to 'genes' if present
if ("GeneName" %in% colnames(dat)) {
  colnames(dat)[colnames(dat) == "GeneName"] <- "genes"
}

# Create a new object to avoid overwriting
dat_subset <- dat

# Reshape to long format
library(dplyr)
library(tidyr)

dat_subset <- dat %>%
  rownames_to_column(var = "genes")

# Step 2: Reshape the data from wide to long format
dat_reshaped <- dat_subset %>%
  pivot_longer(
    cols = -genes,
    names_to = "samples",
    values_to = "TPM"
  )

# Step 3: Join with metadata (ensure column names match)
dat_reshaped <- dat_reshaped %>%
  left_join(metadata.modified, by = c("samples" = "sample"))

# Step 4: View and save output
head(dat_reshaped)
write.csv(dat_reshaped, "dat_reshaped.csv", row.names = FALSE)

# Group and summarize using correct column name (diagnosis instead of treatment)
dat_grouped <- dat_reshaped %>%
  group_by(genes, diagnosis) %>%
  summarise(
    mean_TPM = mean(TPM, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_TPM))  # Or sort by diagnosis if needed

# Save summarized result
write.csv(dat_grouped, "dat_grouped.csv", row.names = FALSE)
