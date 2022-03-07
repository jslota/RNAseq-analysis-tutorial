###MMID coding workshop
#
#RNA seq data analysis
#
#Script 5
#Common data visualizations
#
#2022-03-09
#Jessy Slota

library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(dplyr)

#The volcano plot
#Load differential expression results from terminal timepoint
res <- read.csv("DE results/RML_terminal_DE_results.csv")

#a basic volcano plot
ggplot(res, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point()

#a nicer volcano plot
ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=stat)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey50") +
  geom_vline(xintercept = 0.85, linetype="dashed", color="grey50") +
  geom_vline(xintercept = -0.85, linetype="dashed", color="grey50") +
  scale_color_gradient2(low = "navy", high = "firebrick", mid="grey95", midpoint = 0) +
  theme_classic()

#Plotting number of DE genes at each timepoint
res <- data.frame(timepoint=factor(c("4_wpi","4_wpi","8_wpi","8_wpi","12_wpi","12_wpi","14_wpi","14_wpi","16_wpi","16_wpi","18_wpi","18_wpi","20_wpi","20_wpi","terminal","terminal"),
                                   levels = c("4_wpi","8_wpi","12_wpi","14_wpi","16_wpi","18_wpi","20_wpi","terminal")),
                  direction=rep(c("up", "down"), 8),
                  ngenes=NA)
for (i in c("4_wpi","8_wpi","12_wpi","14_wpi","16_wpi","18_wpi","20_wpi","terminal")) {
  tmp <- read.csv(paste0("DE results/RML_", i, "_DE_results.csv"))
  res[res$timepoint==i&res$direction=="up",]$ngenes <- tmp %>% filter(padj < 0.05, log2FoldChange > 0.85, baseMean > 15) %>% NROW()
  res[res$timepoint==i&res$direction=="down",]$ngenes <- tmp %>% filter(padj < 0.05, log2FoldChange < -0.85, baseMean > 15) %>% NROW()
  rm(tmp)
}

#a basic plot
ggplot(res, aes(x=timepoint, y=ngenes, color=direction, shape=direction)) +
  geom_point()

#a nicer plot
ggplot(res, aes(x=timepoint, y=ngenes, color=direction, shape=direction, label=ngenes)) +
  geom_point() +
  geom_text(nudge_y = 50) +
  scale_color_manual(values=c("navy", "firebrick")) +
  theme_classic()

#The heatmap
#We will use DE genes at terminal timepoint
genes <- read.csv("DE results/RML_terminal_DE_results.csv") %>%  filter(padj < 0.05, abs(log2FoldChange) > 0.85, baseMean > 15) %>% pull(X)

#We will use normalized read-counts to calculate z-scores
zscores <- read.csv("raw data/normalized_read_counts.csv", row.names = 1)
zscores <- as.matrix(zscores[genes,])
zscores <- (zscores-rowMeans(zscores))/matrixStats::rowSds(zscores)

#basic heatmap with hierarchical clustering
pheatmap(zscores)

#nicer heatmap
#specify additional variables required by pheatmap
plot_colors <- rev(colorRampPalette(brewer.pal(11,"PuOr"))(100))#colors for mapping to z-scores
column_annotation <- read.csv("sample_info.csv", row.names = 2)[,-1]#annotation for samples
cls <- brewer.pal(8, "Dark2")#colors for qualitative categorization of samples
annotation_colors <- list(`treatment`=c(`RML`="firebrick", `Mock`="navy"),#this list sets the colors for the annotation
                          `timepoint`=c(`4_wpi`=cls[1],`8_wpi`=cls[2],`12_wpi`=cls[3],`14_wpi`=cls[4],
                                        `16_wpi`=cls[5],`18_wpi`=cls[6],`20_wpi`=cls[7],`terminal`=cls[8]))

pheatmap(zscores, color = plot_colors, annotation_col = column_annotation, annotation_colors = annotation_colors,
         show_rownames = FALSE, show_colnames = FALSE, treeheight_row = 25, treeheight_col = 25)


#Plotting the enrichment results
erch <- read.csv("Enrichr results/full_enrichment_results.csv") #load full enrichment results
res <- erch %>% #filter to top 10 enriched WikiPathways increased at terminal timepoint
  filter(timepoint=="terminal", direction=="up", database == "WikiPathway_2021_Human") %>% 
  arrange(Adjusted.P.value) %>% 
  dplyr::slice(1:10)

#basic enrichment plot
ggplot(res, aes(x=-log10(Adjusted.P.value), y=Term)) +
  geom_point()

#nicer plot
#order Terms based on P-value by converting to a factor
res$Term <- factor(res$Term, levels = res$Term)
ggplot(res, aes(x=-log10(Adjusted.P.value), y=Term, color=Combined.Score, label=Overlap)) +
  geom_point(size=3, alpha=0.5) +
  geom_text(nudge_x = 0.5, hjust=0) +
  scale_color_gradientn(colors=brewer.pal(8, "Oranges")[3:8]) +
  theme_classic()

#Next make a plot for decreased genes at terminal timepoint
res <- erch %>% #filter to top 10 enriched cellular components decreased at terminal timepoint
  filter(timepoint=="terminal", direction=="down", database == "GO_Cellular_Component_2021") %>% 
  arrange(Adjusted.P.value) %>% 
  dplyr::slice(1:10) 
res$Term <- factor(res$Term, levels = res$Term)

ggplot(res, aes(x=-log10(Adjusted.P.value), y=Term, color=Combined.Score, label=Overlap)) +
  geom_point(size=3, alpha=0.5) +
  geom_text(nudge_x = 0.5, hjust=0) +
  scale_color_gradientn(colors=brewer.pal(8, "Purples")[3:8]) +
  theme_classic()
