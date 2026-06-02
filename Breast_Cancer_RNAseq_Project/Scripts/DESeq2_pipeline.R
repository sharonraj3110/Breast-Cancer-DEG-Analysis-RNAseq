# Breast Cancer RNA-seq Differential Expression Analysis
# Dataset: GSE183947
# Method: DESeq2 (R)

# ---------------------------
# 1. Load Libraries
# ---------------------------
library(GEOquery)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# ---------------------------
# 2. Load Dataset
# ---------------------------
gse_list <- getGEO("GSE183947", GSEMatrix = TRUE, getGPL = FALSE)
exprSet <- gse_list[[1]]

# ---------------------------
# 3. Extract Data
# ---------------------------
counts <- exprs(exprSet)
metadata <- pData(exprSet)

# ---------------------------
# 4. Create Group Labels
# ---------------------------
metadata$group <- ifelse(
  grepl("tumor", metadata$source_name_ch1, ignore.case = TRUE),
  "tumor",
  "normal"
)

metadata$group <- factor(metadata$group, levels = c("normal", "tumor"))

# ---------------------------
# 5. Build DESeq2 Object
# ---------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ group
)

# ---------------------------
# 6. Filter Low Counts
# ---------------------------
dds <- dds[rowSums(counts(dds)) > 10, ]

# ---------------------------
# 7. Run DE Analysis
# ---------------------------
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "tumor", "normal"))

# ---------------------------
# 8. Save Results
# ---------------------------
resOrdered <- res[order(res$padj), ]
write.csv(as.data.frame(resOrdered), "DEG_results.csv")

# ---------------------------
# 9. Visualization
# ---------------------------

# Volcano Plot
EnhancedVolcano(res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Breast Cancer vs Normal"
)

# Heatmap (Top 30 genes)
topGenes <- head(rownames(resOrdered), 30)
mat <- assay(vst(dds))[topGenes, ]
mat <- mat - rowMeans(mat)

pheatmap(mat,
         annotation_col = metadata["group"],
         scale = "row")

# PCA Plot
vsd <- vst(dds)
plotPCA(vsd, intgroup = "group")