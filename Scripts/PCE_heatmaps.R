setwd("./sRNA-seq")

sRNA_class = commandArgs(trailingOnly = T)

library(xlsx)

df = "PCE_Table S2 - Expression of sRNAs.xlsx"

if (sRNA_class == "miRNA"){
  sindex = 4
  nclust = 9
  filename = "Temperature_regulated_miRNAs.png"
  height = 7
} else if (sRNA_class == "phasiRNA") {
  sindex = 5
  nclust = 5
  filename = "Temperature_regulated_phasiRNAs.png"
  height = 5
} else if (sRNA_class == "hcsiRNA") {
  sindex = 6
  nclust = 4
  filename = "Temperature_regulated_hcsiRNAs.png"
  height = 6
}

Sample.names = function() {
  tissues = c("Seedling", "Root", "Leaf", "Flower")
  conditions = c("15", "21", "27")
  unit = "°C"
  labels = c()
  for(tissue in tissues) {
    for(cond in conditions) {
      label = paste(tissue, cond, unit, sep = " ")
      labels = c(labels, label)
    }
  }
  return(labels)
}

a = read.xlsx2(df, sheetIndex = sindex, startRow = 3, colIndex = c(1,3:14), header = F, colClasses = c("character", rep("numeric", 12)))
rownames(a) = a[,1]
clabels = Sample.names()
a = a[2:nrow(a), 2:13]
a = as.matrix(a)

library(pheatmap)
library(RColorBrewer)

png(filename, width = 6, height = height, units="in", res=600, pointsize=10, type="cairo-png")

pheatmap(log2(a + 0.01), scale = "row",
         labels_col = clabels,
         fontsize = 8,
         cluster_cols = F, cluster_rows = T,
         cutree_rows = nclust,
         cellwidth = 20, cellheight = 8,
         color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)),
         gaps_col = c(3, 6, 9)
         )
dev.off()
rm(list = ls())