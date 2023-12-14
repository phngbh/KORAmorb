#!/usr/bin/env Rscript

## USAGE: Rscript processing_methylation.R <working dir> <phenotype file> <sample folder> <data file>

args = commandArgs(trailingOnly=TRUE)
setwd(as.character(args[[1]]))

libs <- c("/lustre/groups/cbm01/workspace/phong.nguyen/KORAmorb/Methylation/Rlib", .libPaths())
.libPaths(libs)

library(ggplot2, quietly = TRUE, lib.loc = "Rlib")
library(dplyr, quietly = TRUE, lib.loc = "Rlib")
library(haven, quietly = TRUE, lib.loc = "Rlib")
library(matrixStats, lib.loc = "Rlib")
library(MatrixGenerics, lib.loc = "Rlib")
library(SummarizedExperiment, lib.loc = "Rlib")
library(bumphunter, lib.loc = "Rlib")
library(minfi, lib.loc = "Rlib")
library(DMRcate, lib.loc ="Rlib")
library(ChAMPdata, quietly = TRUE, lib.loc = "Rlib")
library(ChAMP, quietly = TRUE, lib.loc = "Rlib")
library(limma, quietly = TRUE, lib.loc = "Rlib")
library(ebGSEA, lib.loc = "Rlib")
library(fgsea, lib.loc = "Rlib")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19, lib.loc = "Rlib")
library(tidyverse, lib.loc = "Rlib")
library(caret, lib.loc = "Rlib")
library(glmnet, lib.loc = "Rlib")
library(pROC, lib.loc = "Rlib")
#library(globaltest)
library(reactome.db, lib.loc = "Rlib")
library(org.Hs.eg.db, lib.loc = "Rlib")
library(methylGSA, lib.loc = "Rlib")

## Load technical variables
#loadRData <- function(fileName){
#  #loads an RData file, and returns it
#  load(fileName)
#  get(ls()[ls() != "fileName"])
#}
##mvalue = readRDS("mvalue.rds")
#batch = loadRData("./Technical_variables/KF4_PlateChip_1727.RData")
#pcs = loadRData("./Technical_variables/control_probe_pcs_n1727.RData")
#pcs.df = as.data.frame(pcs[,1:20]) %>% mutate(ID = rownames(.))
#batch$ZZ_nr = as.character(batch$ZZ_nr)
#batch.pcs = inner_join(batch, pcs.df, by = c("ZZ_nr" = "ID"))
#cell.type = read.csv("./Houseman/KF4_QN_BMIQ_estimated_cell_distribution_meanimpute473_lessThanOneFALSE.csv", sep = ";")
#cell.type$ZZ_nr = as.character(cell.type$ZZ_nr)
#info = left_join(batch.pcs, cell.type) %>% column_to_rownames("ZZ_nr")
#info$Chip = as.factor(info$Chip)
#info$Batch = as.factor(info$Batch)
#info = rownames_to_column(info, "Sample_Name")

cat("Load necessary files and prepare for statistical analysis\n")
cat("...Load phenotype file\n")
# Load phenotype file
phenotype <- readRDS(as.character(args[[2]]))
target <- colnames(phenotype)[1]
cat("......Phenotype file loaded with dimensions: ",dim(phenotype),"\n")

cat("...Load pathways\n")
# Load gene set data
#lowest_level_pathways = readRDS("lowest_level_pathways.rds")
data("dualmap450kEID")
#geneset_reactome = reactomePathways(names(mapEIDto450k.lv))
#geneset_reactome <- geneset_reactome[intersect(names(geneset_reactome), lowest_level_pathways)]
pathway_list <- readRDS("pathway_list_meth.rds")

cat("...Load data file\n")
# Load processed (filtered and batch corrected) data
beta_processed = readRDS(as.character(args[[4]]))
cat("......Data table loaded with dimensions: ", dim(beta_processed),"\n")

cat("...Load samples\n")
# Load sample file
fs_samples <- readRDS(paste0(args[[3]],"/samples_fs.rds"))
cat("......There are: ", length(fs_samples), " samples loaded\n")

cat("...Extract data and make model matrix for feature selection\n")
# Extract data and model matrix for DMA
int <- intersect(fs_samples, intersect(colnames(beta_processed),rownames(phenotype)))
beta_fs <- beta_processed[,int]
cat("......There are: ", ncol(beta_fs), " samples in whole feature selection set\n")
phenotype <- phenotype[int,,drop=FALSE]
cat("......Dimensions of filtered phenotype table: ",dim(phenotype),"\n")
phenotype[,target] <- as.factor(phenotype[,target])
f <- reformulate(termlabels = target, intercept = F)
modmatrix <- model.matrix(f, data = phenotype)
# modmatrix <- model.matrix(~ 0 + ., data = phenotype)
cat("......Dimensions of filtered model matrix: ",dim(modmatrix),"\n")

# Initiate result list
res = list(edge = list(), ilmn = list(), roc = list(), auc = list(), pathway = list())

# Get subset number
i <- as.numeric(args[[5]])

cat("Start feature selection for resample ",i,"\n")

cat("...Extract relevant samples\n")
subset <- readRDS(paste0(args[[3]],"/samples_",i,".rds"))
int <- intersect(subset, colnames(beta_fs))
beta_sub = beta_fs[,int]
modmatrix_sub <- modmatrix[colnames(beta_sub),,drop = FALSE]
cat("......Subset has ", length(int), " samples\n")

# Differential methylation analysis
cat("...DM analysis\n")
if( any(colnames(beta_sub) != rownames(modmatrix_sub))){
  stop("......columns of data do not match rows of design matrix => Stop!")
}
if( any(table(modmatrix_sub[,1]) < 2) | any(table(modmatrix_sub[,2]) < 2)){
  stop("......One of the factor levels has less than 2 observations => Stop!\n")
}
fit = lmFit(log(beta_sub), modmatrix_sub)
contrast = makeContrasts(paste0(target,"1-",target,"0"), levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contrast)
tmp <- eBayes(tmp)
topdm <- topTable(tmp, sort.by = "P", n = Inf)
sigdm <- rownames(topdm[topdm$P.Value < 2.4e-07,])
cat("......There are ", length(sigdm), " significant differentially methylated probes\n")
saveRDS(sigdm,paste0("fsRes/sigDM_",target,"_",i,".rds"))
ranklist = topdm$P.Value
names(ranklist) = rownames(topdm)
ranklist = sort(ranklist)

# Gene set enrichment analysis
cat("...GSEA\n")
gseaRes <- methylRRA(cpg.pval = ranklist, method = "GSEA", GS.list = pathway_list, GS.idtype = "ENTREZID", minsize = 3, maxsize = 200) %>%
  dplyr::filter(padj < 0.05)
if (nrow(gseaRes) == 0){
  cat("......There is no significantly enriched pathway => Stop!\n")
  res[["edge"]] = NULL
  res[["ilmn"]] = NULL
  res[["pathway"]] = NULL
} else {
  edge_tmp <- gseaRes$core_enrichment %>% paste(collapse = "/") %>% strsplit(split = "/") %>% unlist() %>% unique()
  eid <- mapIdsList(x=org.Hs.eg.db, keys=edge_tmp,keytype="SYMBOL", column="ENTREZID") %>% unlist() %>% unname()
  ilmn <- mapEIDto450k.lv[eid] %>% unlist() %>% unique()
  ilmn <- intersect(ilmn, rownames(beta_processed))
  res[["edge"]] = edge_tmp
  res[["ilmn"]] = ilmn
  res[["pathway"]] = gseaRes$ID
  
  cat("...Lasso\n")
  probelist_tmp = res$ilmn
  oob = setdiff(colnames(beta_fs), colnames(beta_sub))
  x_train = beta_sub[probelist_tmp,] %>% t()
  y_train = phenotype[colnames(beta_sub),target] %>% droplevels()
  up = upSample(x = x_train, y = y_train)
  x_train_up = up[,-ncol(up)]
  y_train_up = up$Class
  x_test = beta_fs[probelist_tmp,oob] %>% t()
  y_test = phenotype[oob,target] %>% droplevels()
  if (nlevels(y_train_up) < 0 | nlevels(y_test) < 2) {
    gsea_inc1$roc = NA
    gsea_inc1$auc = NA
    cat("......The labels in train or test set has only one level => Stop!\n")
  } else {
    fit_tmp = cv.glmnet(x = as.matrix(x_train_up), y = y_train_up, alpha = 1, family ="binomial",
                        nfolds = 3, type.measure = "auc")
    pred_prob = predict(fit_tmp, x_test, s = "lambda.min", type = "response")
    roc <- roc(response = y_test, predictor = pred_prob[,1], levels = c("0","1"))
    auc = auc(roc)
    res[["roc"]] = roc
    res[["auc"]] = auc
  }
}

saveRDS(res,paste0("fsRes/fsRes_",target,"_",i,".rds"))


