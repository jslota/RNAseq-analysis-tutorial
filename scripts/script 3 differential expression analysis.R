###MMID coding workshop
#
#RNA seq data analysis
#
#Script 3
#Differential expression analysis
#
#2022-03-09
#Jessy Slota

library(DESeq2)

#load raw data
read_counts <- read.csv("raw data/raw_read_counts.csv", row.names = 1)#load read count files
sample_info <- read.csv("sample_info.csv", row.names = 2, stringsAsFactors = TRUE)[,-1]#load sample info and set rownames to sample name

#Only keep samples from terminal timepoint
samples <- rownames(sample_info[sample_info$timepoint=="terminal",])

#make DEseq data object
dds <- DESeqDataSetFromMatrix(countData = read_counts[,samples], colData = sample_info[samples,], design = ~treatment)
dds <- DESeq(dds)

#get differential expression results
resultsNames(dds)
res <- results(object = dds, contrast = c("treatment", "RML", "Mock"))#Contrast = RML vs Mock samples

#clean up results file
res <- res[order(res$padj),]
res <- na.omit(res)
res <-as.data.frame(res)

summary(res$padj < 0.05)#Get summary of statistical significance

#Save differential expression results
if (dir.exists("DE results")==FALSE) { dir.create("DE results") }
write.csv(res, "DE results/RML_terminal_DE_results.csv")

#Advanced - For loop that tests every comparison
for (i in unique(sample_info$timepoint)) {
  #get samples
  samples <- rownames(sample_info[sample_info$timepoint==i,])
  #make DEseq data object
  dds <- DESeqDataSetFromMatrix(countData = read_counts[,samples], colData = sample_info[samples,], design = ~treatment)
  dds <- DESeq(dds)
  
  #get differential expression results
  resultsNames(dds)
  res <- results(object = dds, contrast = c("treatment", "RML", "Mock"))
  
  #clean up results file
  res <- res[order(res$padj),]
  res <- na.omit(res)
  res <-as.data.frame(res)
  write.csv(res, paste0("DE results/RML_", i, "_DE_results.csv"))
  print(paste0("saving file... ", "DE results/RML_", i, "_DE_results.csv"))
}
