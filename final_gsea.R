#!/usr/bin/env Rscript

library(dplyr, lib.loc = "/lustre/groups/cbm01/workspace/phong.nguyen/DSPN/Methylation/Rlib")
library(RobustRankAggreg, lib.loc = "/lustre/groups/cbm01/workspace/phong.nguyen/DSPN/Methylation/Rlib")
library(fgsea, lib.loc = "/lustre/groups/cbm01/workspace/phong.nguyen/DSPN/Methylation/Rlib")

args = commandArgs(trailingOnly=TRUE)

setwd(as.character(args[[2]]))

cat("Extract significant pathways\n")
pwlist <- list()
for (i in 1:100){
  res <- read.table(paste0("magma_geneset/",as.character(args[[1]]),"_",i,".gsa.out"), header = T, sep = "", skip = 4) 
  res$adjP <- p.adjust(res$P, method = "fdr")
  res <- arrange(res, adjP) %>% filter(adjP < 0.1)
  pwlist[[i]] <- as.character(res$FULL_NAME) 
}

pwlist <- pwlist[lapply(pwlist,length)>0]
if (length(pwlist) > 3){
  cat("Aggregate pathways into a single list\n")
  set.seed(993)
  pwlist_agg = aggregateRanks(pwlist)
  pwlist_agg$adjP = pwlist_agg$Score*length(pwlist)
  pwlist_agg$adjP = p.adjust(pwlist_agg$adjP, method = "fdr")
  toppw = filter(pwlist_agg, adjP < 0.05)$Name %>% as.character()
} else {
  cat("There are less than three significant pathways to aggregate\n")
}

if (length(toppw) > 0){
  write.table(data.frame(V1 = toppw), paste0(as.character(args[[1]]),"top_pathways.txt"), col.names = F, row.names = F, quote = F)
  cat("Do final GSEA on most significant pathways\n")
  gene_df <- read.table(paste0("magma_gene/",as.character(args[[1]]),".genes.out"), header = T, sep = "") %>%
    arrange(ZSTAT) 
  ranklist <- gene_df$ZSTAT
  names(ranklist) <- gene_df$GENE
  geneset <- readRDS("geneset.rds")
  geneset <- geneset[toppw]
  set.seed(993)
  fgseaRes <- fgsea(pathways = geneset,
                    stats    = ranklist,
                    #minSize  = 15,
                    maxSize  = 200,
                    nperm = 5000) %>% arrange(pval)
  edge_genes <- fgseaRes$leadingEdge %>% unlist() %>% unique()
  
  if (length(edge_genes) > 0){
    cat("Save leading edge genes and their associated SNPs\n")
    write.table(data.frame(V1 = edge_genes), paste0(as.character(args[[1]]),"_leadingEdge_genes.txt"), col.names = F, row.names = F, quote = F)
    gene_snps <- read.table("gene_annot.txt") %>% filter(V1 %in% edge_genes)
    write.table(data.frame(V1 = unique(gene_snps$V2)), paste0(as.character(args[[1]]),"_leadingEdge_SNPs.txt"), col.names = F, row.names = F, quote = F)
  } else {
    cat("There is no leading edge genes\n")
  }
} else {
  cat("There is no significant pathways after aggregation\n")
}

