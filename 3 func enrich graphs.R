#Buble charts functional enrichment

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

pat.genes <- read.csv2("paternal genes heatmap.csv", header = TRUE)
rownames(pat.genes) <- pat.genes$symbol
paternal.genes.enrich <- pat.genes[,-c(1:3)]

#########

# Generate annotations for rows and columns
annotation_col = data.frame(Cluster = factor(pat.genes$cluster.number))
rownames(annotation_col) = pat.genes$symbol

ann_colors = list(
  Cluster = c("0"="#cdcdcd", "1"="#e98383","2"="#eee2b0" , "3"="#91bc75" ,"4"="#a0deb9" ,"5"="#6ca5b3" ,"6"="#a4a5e2")
)

#plot heatmap
pheatmap(t(paternal.genes.enrich2), fontsize = 7, cluster_rows = TRUE,
         color = colorRampPalette(c("#edecec", "#886a9c"))(50),
         annotation_col = annotation_col,
         cluster_cols = TRUE,
         annotation_colors = ann_colors,legend = FALSE,
         border_color = "black",
         cutree_rows = 2, cutree_cols = 3)
         