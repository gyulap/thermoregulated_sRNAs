setwd("./sRNA-seq/ShortStack_results/DESeq2")

timestamp = format(Sys.time(), "%Y.%m.%d_%H.%M.%S")

filtdir = paste0(getwd(), "/Filtered_", timestamp)

if (!dir.exists(filtdir)) {
  dir.create(filtdir)
  }

fc = 1.5
bm = 10
fdr = 0.05

tissues = c("Seedling", "Root", "Leaf", "Flower")

MA_plot = function(a, fc, bm) {
  upreg = as.character(nrow(a[a$passed == T & sign(a$log2FoldChange) == 1,]))
  downreg = as.character(nrow(a[a$passed == T & sign(a$log2FoldChange) == -1,]))
  plot(log2(a$baseMean), a$log2FoldChange,
       type = "p", pch = ifelse(a$passed == T, "×", "."), cex = 1.5,
       col = ifelse(a$passed == T, rgb(0.7,0.2,0.2, alpha=1.0), rgb(0.275,0.51,0.706, alpha=1.0)),
       xlim = c(-3, 24), ylim = c(-8, 8),
       ann = F, axes = F
       )

  axis(1, at = c(0, 8, 16, 24), labels = c(0, 8, 16, 24))
  axis(2, at = c(-8, -4, 0, 4, 8), labels = c(-8, -4, 0, 4, 8))
  title(xlab = "log2(baseMean)", ylab = "log2(FoldChange)")
  abline(h = 0, col = "black")
  abline(h = log2(fc), col = "red", lty = 2)
  abline(h = -log2(fc), col = "red", lty = 2)
  abline(v = log2(bm), col = "red", lty = 2)
  text(20, 6, labels = bquote("" %up% .(upreg)))
  text(20, -6, labels = bquote("" %down% .(downreg)))
  box(which = "plot", lty = "solid", col = "black")
}

pngfile = paste0(filtdir, "/MA-plot_FC", fc, "_bm", bm, "_fdr", fdr, ".png")
png(pngfile, width=9, height=4, units="in", res=600, pointsize=10, type="cairo-png")
mat = matrix(c(0:5, 7:10, 6, 11:14), nrow = 3, ncol = 5, byrow = T)
layout(mat, widths = rep(c(3.5, 5, 5, 5, 5), 3), heights = c(rep(c(1, 5, 5), 5)))
par(oma = c(0, 0, 0, 1))

for(i in 1:6){
  if(i %in% 1:4){
    txt = tissues[i]
    par(mar = c(0, 4, 0, 0) + 0.1)
  }
  else if(i == 5){
    txt = bquote("15 °C" %->% "21 °C")
    par(mar = c(4, 0, 0, 0) + 0.1)
  }
  else if(i == 6){
    txt = bquote("21 °C" %->% "27 °C")
    par(mar = c(4, 0, 0, 0) + 0.1)
  }
  
  plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", axes = F)
  text(x = 0.5, y = 0.5, labels = txt, cex = 1.5)
}

infiles = list.files(getwd(), recursive=F, pattern = "^DESeq2_results_.*[.]txt$")
infiles = infiles[c(4,3,2,1,8,7,6,5)]

par(mar = c(4, 4, 0, 0) + 0.1, mgp = c(2, 1, 0), las = 1)

for(i in 1:8){
  a = read.table(paste0(getwd(), "/", infiles[i]), sep="\t", header = T, row.names = NULL, fill=T, check.names = F, quote = "")
    
  a$passed = ifelse(a$baseMean >= bm
                    & !is.na(a$log2FoldChange)
                    & !is.na(a$padj)
                    & a$padj < fdr
                    & abs(a$log2FoldChange) >= log2(fc),
                    T, F)
  a = a[order(a$passed),]
    
  MA_plot(a = a, fc = fc, bm = bm)
    
  tablefile = paste0(filtdir, "/", gsub(".txt", "", infiles[i]), "_FC", fc, "_bm", bm, "_fdr", fdr, ".txt")
  a = a[order(a$log2FoldChange, decreasing = T),]
  write.table(a[a$passed == T,], file=tablefile, sep = "\t", quote = F, row.names = F, col.names = T)
  }

dev.off()
rm(list=ls())