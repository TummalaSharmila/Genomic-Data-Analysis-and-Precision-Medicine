# Install required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db"))
install.packages(c("dplyr", "writexl"))

# Load required libraries
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(writexl)

# Set working directory
setwd("/Users/sharmilatummala/Documents/Precision MEdicine")

# ---- 1. Read in data ----
counts <- read.csv("/Users/sharmilatummala/Downloads/gene_counts.csv", row.names = 1, header = TRUE)
meta_full <- read.csv("/Users/sharmilatummala/Downloads/assignment_2_info.csv", header = TRUE)

# ---- 2. Clean and extract metadata ----
meta <- meta_full[, c("Condition", "Time.Point", "SRA.Accession")]
colnames(meta) <- c("Condition", "Time_Point", "SRA_Accession")
meta <- meta[match(colnames(counts), meta$SRA_Accession), ]

# ---- 3. Convert grouping variables to factors and clean names ----
meta$Condition <- as.factor(meta$Condition)
meta$Time_Point <- as.factor(meta$Time_Point)
levels(meta$Condition) <- make.names(levels(meta$Condition))
levels(meta$Time_Point) <- paste0("X", levels(meta$Time_Point))
rownames(meta) <- meta$SRA_Accession

# ---- 4. DESeq2 analysis: all samples, Condition + Time_Point ----
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ Condition + Time_Point)
dds <- DESeq(dds)

# ---- 5. DEGs: Mock vs SARS-CoV-2 (all time points) ----
res_mock_vs_SARS <- results(dds, contrast = c("Condition", "Mock.", "SARS.CoV.2"))
sig_mock_vs_SARS <- subset(res_mock_vs_SARS, padj < 0.05)
write.csv(sig_mock_vs_SARS, "DEGs_Mock_vs_SARS.csv")

# ---- 6. DEGs: SARS-CoV-2 at 24H vs SARS-CoV-2 at 72H ----
meta_SARS <- meta[meta$Condition == "SARS.CoV.2", ]
counts_SARS <- counts[, meta_SARS$SRA_Accession]
dds_SARS <- DESeqDataSetFromMatrix(countData = counts_SARS, colData = meta_SARS, design = ~ Time_Point)
dds_SARS <- DESeq(dds_SARS)
res_24H_vs_72H <- results(dds_SARS, contrast = c("Time_Point", "X24", "X72"))
sig_24H_vs_72H <- subset(res_24H_vs_72H, padj < 0.05)
write.csv(sig_24H_vs_72H, "DEGs_SARS_24H_vs_72H.csv")

# ----  Extract normalized expression matrix for all genes and samples ----
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, "Expression_Levels_All_Genes_Per_Sample_Annotated_With_Conditions.csv")

# ---- 7. GO Term Enrichment Analysis ----
# For Mock vs SARS-CoV-2 DEGs
gene_ids_mock_SARS <- rownames(sig_mock_vs_SARS)
gene_list_mock_SARS <- bitr(gene_ids_mock_SARS, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
go_mock_SARS <- enrichGO(gene = gene_list_mock_SARS$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH")
write.csv(go_mock_SARS, "GO_Enrichment_Mock_vs_SARS.csv")

# For SARS-CoV-2 24H vs 72H DEGs
gene_ids_24H_72H <- rownames(sig_24H_vs_72H)
gene_list_24H_72H <- bitr(gene_ids_24H_72H, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
go_24H_72H <- enrichGO(gene = gene_list_24H_72H$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH")
write.csv(go_24H_72H, "GO_Enrichment_SARS_24H_vs_72H.csv")

cat("Analysis complete. Output files saved.\n")


# Loading required libraries
library(dplyr)

# Converting DESeqResults to data frame
sig_mock_vs_SARS_df <- as.data.frame(sig_mock_vs_SARS)
sig_24H_vs_72H_df <- as.data.frame(sig_24H_vs_72H)

# Summarizing DEG results for Mock vs SARS-CoV-2
summary_mock_vs_SARS <- sig_mock_vs_SARS_df %>%
  summarise(
    Total_DEGs = n(),
    Upregulated = sum(log2FoldChange > 0, na.rm = TRUE),
    Downregulated = sum(log2FoldChange < 0, na.rm = TRUE),
    Mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    Median_log2FC = median(log2FoldChange, na.rm = TRUE),
    Min_padj = min(padj, na.rm = TRUE),
    Max_padj = max(padj, na.rm = TRUE)
  )

# Summarizing DEG results for SARS-CoV-2 24H vs 72H
summary_24H_vs_72H <- sig_24H_vs_72H_df %>%
  summarise(
    Total_DEGs = n(),
    Upregulated = sum(log2FoldChange > 0, na.rm = TRUE),
    Downregulated = sum(log2FoldChange < 0, na.rm = TRUE),
    Mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    Median_log2FC = median(log2FoldChange, na.rm = TRUE),
    Min_padj = min(padj, na.rm = TRUE),
    Max_padj = max(padj, na.rm = TRUE)
  )

# Print summaries
print("Summary: Mock vs SARS-CoV-2")
print(summary_mock_vs_SARS)

print("Summary: SARS-CoV-2 24H vs 72H")
print(summary_24H_vs_72H)

# Optionally, save summaries to CSV
write.csv(summary_mock_vs_SARS, "Summary_Mock_vs_SARS.csv")
write.csv(summary_24H_vs_72H, "Summary_SARS_24H_vs_72H.csv")


