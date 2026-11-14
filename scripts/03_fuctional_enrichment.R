# ==============================================================================
# Script 03: Gene Annotation & GO/KEGG Functional Enrichment
# ==============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# 1. Load DESeq2 results
res_df <- readRDS("results/res_df.rds")

# Strip Ensembl version suffixes if present (e.g., ENSG00000171195.12 -> ENSG00000171195)
rownames(res_df) <- sub("\\..*$", "", rownames(res_df))

# 2. Annotate top significant genes with official HGNC symbols
top_ids <- rownames(head(res_df, 20))
gene_map <- bitr(top_ids, fromType = "ENSEMBL", toType = c("SYMBOL", "GENENAME"), OrgDb = org.Hs.eg.db)
write.csv(gene_map, "results/Top20_Annotated_Genes.csv", row.names = FALSE)

# 3. Extract significant Upregulated & Downregulated gene lists (padj < 0.05)
up_ensembl   <- rownames(subset(res_df, padj < 0.05 & log2FoldChange > 0.5))
down_ensembl <- rownames(subset(res_df, padj < 0.05 & log2FoldChange < -0.5))

up_entrez   <- bitr(up_ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
down_entrez <- bitr(down_ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# 4. Run GO Over-Representation Analysis (BP, CC, MF)
run_go_and_plot <- function(gene_list, ont_type, direction_label, file_ext = "png") {
  ego <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    ont           = ont_type,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  if (!is.null(ego) && nrow(ego) > 0) {
    write.csv(as.data.frame(ego), paste0("results/GO_", ont_type, "_", direction_label, "_GSE349730.csv"), row.names = FALSE)
    p <- dotplot(ego, showCategory = 15) + ggtitle(paste0("GO ", ont_type, " - ", direction_label))
    ggsave(paste0("plots/GO_", ont_type, "_dotplot_", direction_label, "_GSE349730.", file_ext),
           plot = p, width = 9, height = 7, dpi = 300)
  }
}

run_go_and_plot(up_entrez, "BP", "Up", "jpg")
run_go_and_plot(up_entrez, "CC", "Up", "png")
run_go_and_plot(up_entrez, "MF", "Up", "jpg")
run_go_and_plot(down_entrez, "BP", "Down", "png")

# 5. KEGG Pathway Enrichment
ekegg_up <- enrichKEGG(gene = up_entrez, organism = "hsa", pvalueCutoff = 0.05)
if (!is.null(ekegg_up) && nrow(ekegg_up) > 0) {
  write.csv(as.data.frame(ekegg_up), "results/KEGG_Up_GSE349730.csv", row.names = FALSE)
}

cat("Functional enrichment complete!\n")
