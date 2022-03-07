###MMID coding workshop
#
#RNA seq data analysis
#
#Script 2
#Normalizing the data and assessing variation
#
#2022-03-09
#Jessy Slota

library(DESeq2)
library(ggplot2)
library(RColorBrewer)

#load raw data
read_counts <- read.csv("raw data/raw_read_counts.csv", row.names = 1)#load read count files
sample_info <- read.csv("sample_info.csv", row.names = 2, stringsAsFactors = TRUE)[,-1]#load sample info and set rownames to sample name

summary(colnames(read_counts)==rownames(sample_info))#make sure samples are in order

#make DEseq data object
dds <- DESeqDataSetFromMatrix(countData = read_counts, colData = sample_info, design = ~treatment+timepoint)
dds <- DESeq(dds)

#plot dispersion estimates to examine normalization
plotDispEsts(dds)

#Get normalized counts and make PCA plots
norm_counts <- vst(dds)#extract normalized read counts
plotPCA(norm_counts, intgroup=c("treatment", "timepoint"))#make a basic PCA plot

#make a custom PCA plot with ggplot
pcaData <- plotPCA(norm_counts, intgroup=c("treatment", "timepoint"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=timepoint, shape=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = brewer.pal(8, "Dark2")) +
  coord_fixed() +
  theme_classic()

#Save normalized read counts for visualization later
write.csv(assay(norm_counts), "raw data/normalized_read_counts.csv")
