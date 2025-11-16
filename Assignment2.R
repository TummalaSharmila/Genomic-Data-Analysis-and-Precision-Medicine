# Install required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db"))

# Install CRAN packages
install.packages(c("dplyr", "writexl"))
install.packages("writexl")
# Load required libraries
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(writexl)

# Set working directory to the desired path
setwd("/Users/sharmilatummala/Documents/Precision MEdicine")

# ---- 1. Read in data ----

# Read master count matrix and metadata
counts <- read.csv("/Users/sharmilatummala/Downloads/gene_counts.csv", row.names = 1, header = TRUE)        # Update filename if needed
meta_full <- read.csv("/Users/sharmilatummala/Downloads/assignment_2_info.csv", header = TRUE)

# ---- 2. Clean and extract metadata ----

# Only keep Condition, Time Point, SRA Accession columns
meta <- meta_full[, c("Condition", "Time.Point", "SRA.Accession")]
colnames(meta) <- c("Condition", "Time_Point", "SRA_Accession")

# Ensure order matches counts (columns = SRA Accession)
meta <- meta[match(colnames(counts), meta$SRA_Accession), ]

# ---- 3. Convert grouping variables to factors and clean names ----

meta$Condition <- as.factor(meta$Condition)
meta$Time_Point <- as.factor(meta$Time_Point)

# Clean factor level names to be safe
levels(meta$Condition) <- make.names(levels(meta$Condition))
levels(meta$Time_Point) <- paste0("X", levels(meta$Time_Point))

meta

rownames(meta) <- meta$SRA_Accession

all(rownames(meta) %in% colnames(counts))
meta


# ---- 4. DESeq2 analysis: all samples, Condition + Time_Point ----

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ Condition + Time_Point)
dds <- DESeq(dds)

# ---- 5. Differential genes: SARS-CoV-2 vs Mock (all time points) ----

# Check final factor level names for Condition
print(levels(meta$Condition))         # Should be like "Mock." "SARS.CoV.2"
res_SARS_vs_Mock <- results(dds, contrast = c("Condition", "SARS.CoV.2", "Mock."))
sig_SARS_vs_Mock <- subset(res_SARS_vs_Mock, padj < 0.05)
write.csv(sig_SARS_vs_Mock, "Significant_Genes_SARS_vs_Mock.csv")

# ---- 6. Differential genes: 72H vs 24H (SARS-CoV-2 only) ----

meta_SARS <- meta[meta$Condition == "SARS.CoV.2", ]
counts_SARS <- counts[, meta_SARS$SRA_Accession]
dds_SARS <- DESeqDataSetFromMatrix(countData = counts_SARS, colData = meta_SARS, design = ~ Time_Point)
dds_SARS <- DESeq(dds_SARS)
print(levels(meta_SARS$Time_Point))   # Should be "X24", "X72"
res_72H_vs_24H <- results(dds_SARS, contrast = c("Time_Point", "X72", "X24"))
sig_72H_vs_24H <- subset(res_72H_vs_24H, padj < 0.05)
write.csv(sig_72H_vs_24H, "Significant_Genes_72H_vs_24H.csv")

# ---- 7. Extract normalized expression matrix for all genes and samples ----

norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, "Expression_Levels_All_Genes_Per_Sample_Annotated_With_Conditions.csv")
# ---- 8. GO Enrichment analysis ----

# For DEGs: SARS-CoV-2 vs Mock
gene_ids_SARS <- rownames(sig_SARS_vs_Mock)
gene_list_SARS <- bitr(gene_ids_SARS, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
go_SARS_vs_Mock <- enrichGO(gene = gene_list_SARS$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH")
write.csv(go_SARS_vs_Mock, "GO_Enrichment_SARS_vs_Mock.csv")

# For DEGs: SARS-CoV-2 72H vs 24H
gene_ids_72H <- rownames(sig_72H_vs_24H)
gene_list_72H <- bitr(gene_ids_72H, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
go_72H_vs_24H <- enrichGO(gene = gene_list_72H$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH")
write.csv(go_72H_vs_24H, "GO_Enrichment_72H_vs_24H.csv")

cat("Analysis complete. Output files saved.\n")



# Loading required libraries
library(dplyr)

# Converting DESeqResults to data frame
sig_SARS_vs_Mock_df <- as.data.frame(sig_SARS_vs_Mock)
sig_72H_vs_24H_df <- as.data.frame(sig_72H_vs_24H)

# Summarizing DEG results for SARS-CoV-2 vs Mock
summary_SARS_vs_Mock <- sig_SARS_vs_Mock_df %>%
  summarise(
    Total_DEGs = n(),
    Upregulated = sum(log2FoldChange > 0, na.rm = TRUE),
    Downregulated = sum(log2FoldChange < 0, na.rm = TRUE),
    Mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    Median_log2FC = median(log2FoldChange, na.rm = TRUE),
    Min_padj = min(padj, na.rm = TRUE),
    Max_padj = max(padj, na.rm = TRUE)
  )

# Summarizing DEG results for 72H vs 24H (SARS-CoV-2 only)
summary_72H_vs_24H <- sig_72H_vs_24H_df %>%
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
print("Summary: SARS-CoV-2 vs Mock")
print(summary_SARS_vs_Mock)

print("Summary: 72H vs 24H (SARS-CoV-2 only)")
print(summary_72H_vs_24H)

# Optionally, save summaries to CSV
write.csv(summary_SARS_vs_Mock, "Summary_SARS_vs_Mock.csv")
write.csv(summary_72H_vs_24H, "Summary_72H_vs_24H.csv")


