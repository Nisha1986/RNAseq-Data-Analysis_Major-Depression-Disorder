library(ggplot2)
library(tidyverse)
#data
dat.long<-read.csv("dat_grouped.csv")
dim(dat.long)
head(dat.long)
dat.reshaped<-read.csv("dat_reshaped.csv")
# Load the mapping file
gene_map <- read.csv("cleaned_mdd_dataset_with_genesymbols.csv", stringsAsFactors = FALSE)
# Preview the mapping
head(gene_map)
dim(gene_map)
#dat_grouped to be used generated frpm previous demo
# Step 2: Map gene names to dat.reshaped
# Step 2: Join gene mapping to the main data
dat.reshaped_mapped <- dat.reshaped %>%
  left_join(gene_map, by = c("genes" = "ensembl_gene_id"))

# Step 3: Replace 'genes' column with gene name (fallback to Ensembl if missing)
dat.reshaped_mapped <- dat.reshaped_mapped %>%
  mutate(genes = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "", genes, hgnc_symbol)) %>%
  select(-hgnc_symbol)  # Drop hgnc_symbol column as it's now merged into 'genes'

# Optional: Preview the data
head(dat.reshaped_mapped)
library(ggplot2)

ggplot(dat.reshaped_mapped, aes(x = diagnosis, y = TPM, fill = diagnosis)) +
  geom_boxplot() +
  labs(title = "TPM Expression by Diagnosis", x = "Diagnosis", y = "TPM") +
  theme_minimal()

# 1️⃣ Boxplot: TPM distribution across samples
dat_wide <- dat.reshaped_mapped %>%
  select(genes, samples, TPM) %>%
  pivot_wider(names_from = samples, values_from = TPM)
library(ggplot2)
ggplot(dat.reshaped_mapped, aes(x = samples, y = TPM, color = diagnosis)) +
  geom_point(size = 2) +
  labs(title = "TPM Expression per Sample",
       x = "Sample", y = "TPM") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 2️⃣ Density Plot: Distribution of TPM expression per sample
ggplot(dat.reshaped_mapped, aes(x = TPM, fill = Sample)) +
  geom_density(alpha = 0.4) +
  labs(title = "Density Plot: TPM Distribution per Sample",
       x = "samples", y = "TPM") +
  theme_minimal() +
  theme(legend.position = "none")
# Step 4: Save the updated data
write.csv(dat.reshaped_mapped, "dat_reshaped_with_GeneNames.csv", row.names = FALSE)
# Define your top genes (example list)
top_genes <- c("OR7A5", "SCML2P2", "FAT2", "SERPINE1", "PRR35", "RN7SL79P", 
               "UBE2CP2", "SLC9A3-AS1", "GAPDHP1", "MSLN", "C1QL2", "SLC25A41", 
               "ENSG00000228150", "ENSG00000249320", "ENSG00000233052")
# Filter the data for those genes
dat_top <- dat.reshaped_mapped %>% filter(genes %in% top_genes)
# Basic grouped bar plot
ggplot(dat_top, aes(x = samples, y = TPM, fill = diagnosis)) +
  geom_col(position = "dodge") +
  facet_wrap(~ genes, scales = "free_y") +
  theme_bw() + 
  labs(
    title = "TPM Expression of Top Genes Across Samples",
    x = "Sample ID",
    y = "TPM"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, size = 7, hjust = 1),
    legend.position = "top"
  )
ggsave("TPM_TopGenes_Barplot.png", width = 14, height = 8, dpi = 300)

#densityplot
# Add log2 TPM for better visualization
dat_top <- dat_top %>%
  mutate(logTPM = log2(TPM + 1))

# Plot density for each gene
ggplot(dat_top, aes(x = logTPM, fill = diagnosis)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ genes, scales = "free") +
  theme_minimal() +
  labs(
    title = "log2(TPM + 1) Density Plot for Top 15 Genes",
    x = "log2(TPM + 1)",
    y = "Density"
  ) +
  scale_fill_manual(values = c("Normal" = "#4F81BD", "MDD" = "#C0504D")) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)
  )

# Save the plot
ggsave("Top15Genes_TPM_DensityPlot.png", width = 14, height = 8, dpi = 300)

#boxplot
ggplot(dat_top, aes(x = diagnosis, y = TPM, fill = diagnosis)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~ genes, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "TPM Expression of Top Genes (Box Plot)",
    x = "Diagnosis",
    y = "TPM"
  ) +
  scale_fill_manual(values = c("Normal" = "#4F81BD", "MDD" = "#C0504D")) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

ggsave("TopGenes_BoxPlot.png", width = 14, height = 8, dpi = 300)

#scatterplot
# Boxplot
ggplot(dat_top, aes(x = diagnosis, y = TPM, fill = diagnosis)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~ genes, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "TPM Expression of Top Genes (Box Plot)",
    x = "Diagnosis",
    y = "TPM"
  ) +
  scale_fill_manual(values = c("Normal" = "#4F81BD", "MDD" = "#C0504D")) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

ggsave("TopGenes_BoxPlot.png", width = 14, height = 8, dpi = 300)
#Scatter plot
ggplot(dat_top, aes(x = samples, y = TPM, color = diagnosis, group = diagnosis)) +
  geom_point(size = 2) +
  facet_wrap(~ genes, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "TPM Expression of Top Genes (Scatter Plot)",
    x = "Samples",
    y = "TPM"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, size = 7, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  ) +
  scale_color_manual(values = c("Normal" = "#4F81BD", "MDD" = "#C0504D"))

ggsave("TopGenes_ScatterPlot.png", width = 14, height = 8, dpi = 300)
#heatmap
# Pivot to wide format for heatmap: genes as rows, samples as columns
library(pheatmap)
heat_data <- dat_top %>%
  select(genes, samples, TPM) %>%
  pivot_wider(names_from = samples, values_from = TPM)

# Set gene names as row names
heat_mat <- as.data.frame(heat_data)
rownames(heat_mat) <- heat_mat$genes
heat_mat <- heat_mat[ , -1]  # remove the 'genes' column

# Log2 transform
heat_mat_log <- log2(heat_mat + 1)

# Heatmap
pheatmap(heat_mat_log,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Heatmap of Top Genes (log2 TPM)",
         fontsize_row = 8,
         fontsize_col = 7)
ggsave('Topgene_heatmap.png',width = 14, height = 8, dpi = 300)
