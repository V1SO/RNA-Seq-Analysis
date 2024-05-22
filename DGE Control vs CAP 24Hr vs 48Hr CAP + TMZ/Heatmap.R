library(pheatmap)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(viridis)

log2Counts=log2.counts
diffExpressed <- res.sig.TMZ24
countsDiffExpressed=log2Counts[rownames(diffExpressed)[1:100], ]
rownames(countsDiffExpressed)=countsDiffExpressed$symbol
countsDiffExpressed=countsDiffExpressed[, 9:16]
mean24=(rowMeans(countsDiffExpressed[,1:4]))
mean48=(rowMeans(countsDiffExpressed[,5:8]))
countsDiffExpressed=cbind(countsDiffExpressed, mean24)
countsDiffExpressed=cbind(countsDiffExpressed, mean48)
countsDiffExpressed=countsDiffExpressed[,9:10]
#patients=read_csv("../PatientIdentifiers.csv")
#colnames(countsDiffExpressed) <- colnames(patients)
colnames(countsDiffExpressed) <- c("Plasma 24 hour", "Plasma 48 hour + TMZ")

p=pheatmap(countsDiffExpressed)

mat_cluster_cols <- hclust(dist(t(countsDiffExpressed)))
plot(mat_cluster_cols, main = "Unsorted Dendrogram",
     xlab = "", sub = "")
library(dendsort)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

mat_cluster_cols <- sort_hclust(mat_cluster_cols)
plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")
mySampleCol = data.frame(sample=as.factor(c(rep("SLE", 99), rep("CTL", 18))))
rownames(mySampleCol) = colnames(patients)
mat_cluster_rows <- sort_hclust(hclust(dist(countsDiffExpressed)))
newCols <- colorRampPalette(grDev)
my_color = list(sample=c(SLE = "#4B0C6BFF", CTL = "#FB9A06FF"))
png(file="Heatmap502.png", width=2000, height=5000, res=280)
pheatmap(
  mat               = countsDiffExpressed,
  color             = inferno(20),
  cluster_cols      = mat_cluster_cols,
  cluster_rows      = mat_cluster_rows,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  #annotation_colors = mat_colors,
  cellwidth         = 60,
  cellheight        = 10,
  #drop_levels       = TRUE,
  #fontsize          = 14,
  legend_breaks     = 1:20,
  border_color      = NA,
  legend_labels     = c(as.character(1:20)),
  main              = "Sorted Dendrograms",
  cex               = 1,
  annotation_colors = my_color,
  annotion_col      = mySampleCol,
  annotation        = mySampleCol,
  annotation_legend = TRUE,
  annotation_names_col = TRUE,
)
dev.off()
#ggsave("Heatmap2.png", plot = last_plot(), device = "png", dpi = 500, limitsize = FALSE)