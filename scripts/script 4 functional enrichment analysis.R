###MMID coding workshop
#
#RNA seq data analysis
#
#Script 4
#Functional enrichment analysis
#
#2022-03-09
#Jessy Slota

library(enrichR)
library(dplyr)

#Identify DE genes
res <- read.csv("DE results/RML_terminal_DE_results.csv")

#Get some genes to test in enrichr
res %>% filter(padj < 0.05, log2FoldChange > 0.85, baseMean > 15)%>% 
  pull(X) %>% 
  writeClipboard()#copies to clipboard... paste at https://maayanlab.cloud/Enrichr/

#Make list of databases you are interested in
dbs <- c("WikiPathway_2021_Human", "GO_Cellular_Component_2021", "PanglaoDB_Augmented_2021")

#Genes with increased abundance
genes <- res %>% filter(padj < 0.05, log2FoldChange > 0.85, baseMean > 15) %>% pull(X)

#Run Enrichr
enrch <- enrichr(genes, dbs)

#convert results from list to data frame
for (i in dbs) {
  enrch[[i]]$database <- i
}
enrch <- do.call(rbind, enrch)


###Advanced - for loop to get full enrichment results from every comparison
#increased and decreased genes at 8 timepoints = 16 comparisons

#Make list of databases you are interested in
dbs <- c("WikiPathway_2021_Human", "GO_Cellular_Component_2021", "PanglaoDB_Augmented_2021")

#Make empty list for full results
full_enrch_results <- list()
#Loop through every list of genes
sample_info <- read.csv("sample_info.csv", row.names = 2)[,-1]
for (i in unique(sample_info$timepoint)) {
  res <- read.csv(paste0("DE results/RML_", i, "_DE_results.csv"))
  
  ##Genes with increased abundance
  genes <- res %>% filter(padj < 0.05, log2FoldChange > 0.85, baseMean > 15) %>% pull(X)
  if (length(genes) > 0) {
    enrch <- enrichr(genes, dbs)
    for (j in dbs) {
      enrch[[j]]$timepoint <- i
      enrch[[j]]$direction <- "up"
      enrch[[j]]$database <- j
    }
    enrch <- do.call(rbind, enrch)#convert from list to data frame
    full_enrch_results[[paste0(i, "_up")]] <- enrch
    print(paste0("analysis complete... ", i, "_up"))
  }
  
  ##Genes with increased abundance
  genes <- res %>% filter(padj < 0.05, log2FoldChange < -0.85, baseMean > 15) %>% pull(X)
  if (length(genes) > 0) {
        enrch <- enrichr(genes, dbs)
    for (j in dbs) {
      enrch[[j]]$timepoint <- i
      enrch[[j]]$direction <- "down"
      enrch[[j]]$database <- j
    }
    enrch <- do.call(rbind, enrch)#convert from list to data frame
    full_enrch_results[[paste0(i, "_down")]] <- enrch
    print(paste0("analysis complete... ", i, "_down"))
  }
}
full_enrch_results <- do.call(rbind, full_enrch_results)
row.names(full_enrch_results) <- seq(1:nrow(full_enrch_results))

#save final results for later
if (dir.exists("Enrichr results")==FALSE) { dir.create("Enrichr results") }
write.csv(full_enrch_results, "Enrichr results/full_enrichment_results.csv")
