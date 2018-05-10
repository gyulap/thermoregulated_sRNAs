library(dplyr)

setwd("./sRNA-seq/ShortStack_results")

outdir = paste0(getwd(), "/DESeq2")
filtdir = paste0(outdir, "/Filtered")
dir.create(outdir)
dir.create(filtdir)

tissues = c("Seedling", "Root", "Leaf", "Flower")
conditions = c("15", "21", "27")

fc = 1.5
bm = 10
fdr = 0.05

input_data = read.table("Counts.txt", header=T, sep="\t", row.names=2)

library(DESeq2)

normcounts = data.frame("Name" = rownames(input_data), stringsAsFactors = F)
mapping = data.frame("condition" = "", "tissue" = "")
sf = data.frame("Sample" = "", "Factor" = 0)

DESeq2.test = function(ref, cond) {
  results = results(dds, contrast = c("condition", cond, ref),
                    alpha = fdr, lfcThreshold = log2(fc), altHypothesis = "greaterAbs", independentFiltering = T)
  results = lfcShrink(dds, contrast = c("condition", cond, ref), res = results)
  outfilename = paste0(outdir, "/DESeq2_results_", ref, "_to_", cond, "_", tissue, "_", timestamp, ".txt")
  results = as.data.frame(results)
  results$Name = rownames(results)
  results = results[, c(ncol(results), 1:(ncol(results)-1))]
  write.table(results, file = outfilename, sep = "\t", quote = F, row.names = F, col.names = T)
}

for (tissue in tissues) {
  countData = select(input_data, contains(tissue))
  colData = data.frame(
    "tissue" = rep(tissue, 6),
    "condition" = rep(conditions, each = 2),
    row.names = colnames(countData)
  )
  
  mapping = rbind(mapping, colData)
  
  dds = DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
  dds = DESeq(dds)
  
  normcount = as.data.frame(counts(dds, normalized = T))
  normcount = mutate(normcount, "Name" = as.character(rownames(normcount)))
  normcount = select(normcount, "Name", everything())
  normcounts = full_join(normcounts, normcount, by="Name")
  
  factors = data.frame("Factor" = 1/sizeFactors(dds))
  factors$Sample = rownames(factors)
  factors = factors[, c(2, 1)]
  
  sf = rbind(sf, factors)
  
  #Calculating and outputting the results of DE between 15 and 21 °C

  DESeq2.test(ref = "15", cond = "21")
  
  #Calculating and outputting the results of DE between 21 and 27 °C
  
  DESeq2.test(ref = "21", cond = "27")
  
  }

write.table(normcounts,
            file=paste0(outdir, "/DESeq2_norm_counts_", timestamp, ".txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)

write.table(mapping[2:nrow(mapping),],
            file=paste0(outdir, "/DESeq2_norm_colData_", timestamp, ".txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)

write.table(sf[2:nrow(sf),],
            file=paste0(outdir, "/DESeq2_norm_factors_", timestamp, ".txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=ls())