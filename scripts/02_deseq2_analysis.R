# ==============================================================================
# Script 02: Multifactorial DESeq2 Analysis & Core Visualizations
# ==============================================================================

library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)

# 1. Load preprocessed data
filtered_counts <- readRDS("results/filtered_counts.rds")
coldata <- readRDS("results/coldata.rds")

# 2. Construct DESeqDataSet with multifactorial design (~ Patient + Matrix)
dds <- DESeqDataSetFromMatrix(
  countData = round(filtered_counts),
  colData   = coldata,
  design    = ~ Patient + Matrix
)

# 3. Run DESeq2 differential expression
dds <- DESeq(dds)
res <- results(dds, contrast = c("Matrix", "MGCOL", "MG"), alpha = 0.05)
res_ordered <- res[order(res$padj), ]

# Export differential expression tables
res_df <- as.data.frame(res_ordered)
write.csv(res_df, "results/GSE349730_DESeq_all_results.csv")
write.csv(subset(res_df, padj < 0.05 & log2FoldChange > 1), "results/GSE349730_Up_adj.csv")
write.csv(subset(res_df, padj < 0.05 & log2FoldChange < -1), "results/GSE349730_Down_adj.csv")

# 4. Variance-Stabilizing Transformation (VST) & PCA Plot
vsd <- vst(dds, blind = FALSE)
pca_plot <- plotPCA(vsd, intgroup = "Matrix") +
  theme_bw(base_size = 14) +
  ggtitle("PCA - GSE349730 (Matrix Effect)")
ggsave("plots/PCA_GSE349730.png", plot = pca_plot, width = 8, height = 6, dpi = 300)

# 5. EnhancedVolcano Plot
volcano_plot <- EnhancedVolcano(
  res_df,
  lab = rownames(res_df),
  x = "log2FoldChange",
  y = "padj",
  title = "GSE349730: MGCOL vs MG",
  pCutoff = 0.05,
  FCcutoff = 1.0
)
ggsave("plots/Volcano_GSE349730.jpg", plot = volcano_plot, width = 9, height = 7, dpi = 300)

# 6. Top 30 DEGs Heatmap
top_genes <- rownames(head(res_ordered, 30))
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
anno_col <- coldata[, c("Patient", "Matrix")]

png("plots/Heatmap_GSE349730.png", width = 2400, height = 2200, res = 300)
pheatmap(
  mat,
  annotation_col = anno_col,
  show_colnames = TRUE,
  fontsize_col = 7,
  main = "Top DE Genes - GSE349730"
)
dev.off()

saveRDS(res_df, "results/res_df.rds")
cat("DESeq2 analysis and visualizations complete!\n")
