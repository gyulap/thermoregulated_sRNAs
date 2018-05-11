setwd("./RNA-seq/kallisto_results")

library(sleuth)
base_dir = "./kallisto_files"
#ann = read.table("/home/gyulap/Dokumentumok/Arabidopsis/Gyuri_temperature/kallisto/TAIR10_functional_descriptions_20140331.txt",
#                 sep="\t", header = F, fill = T, row.names = NULL, quote = "", stringsAsFactors = F)
sample_id = dir(file.path(base_dir))
kal_dirs = sapply(sample_id, function(id) file.path(base_dir, id))
tissue = rep(c("Flower", "Leaf", "Root", "Seedling"), 1, each = 3)
temperature = rep(c("15", "21", "27"), 4, each = 1)

s2c = data.frame("sample" = sample_id,
                 "tissue" = tissue,
                 "temperature" = temperature,
                 "group" = paste(tissue, temperature, sep = "_"),
                 "path" = kal_dirs, stringsAsFactors = F)

so = sleuth_prep(s2c, ~ group,
                 #target_mapping = ann,
                 extra_bootstrap_summary = T, read_bootstrap_tpm = T)

stat = summary(so)
stat = stat[, c(1,3,2)]
colnames(stat) = c("", "Passed", "Mapped")
stat$Unmapped = stat$Passed - stat$Mapped
write.table(stat, file="kallisto_stat.txt", sep="\t", row.names = F, quote = F)

library(reshape2)
kt = kallisto_table(so, normalized = T, use_filtered = F)
kt2 = dcast(kt, target_id ~ sample, value.var = "tpm")
kt2 = kt2[, c(1,11,12,13,8,9,10,5,6,7,2,3,4)]
write.table(kt2, file="kallisto_isoform_table.txt", sep="\t", row.names = F, quote = F)
rm(list=ls())