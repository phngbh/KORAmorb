#!/usr/bin/env Rscript

library(dplyr)

args = commandArgs(trailingOnly=TRUE)
phenotype <- args[[1]]

perf_validate_list <- list()
perf_test_list <- list()
for (i in 1:100){
  file <- paste0("res_",phenotype,"/",phenotype,"_FFS_glmnet_AUPRC_union_",i,"_fromClinical.rds")
  if(file.exists(file)){
    res <- readRDS(file)
    perf_validate <- res$perf_validate
    perf_validate_list[[i]] <- perf_validate
    perf_test <- res$perf_test
    perf_test_list[[i]] <- perf_test
  } else {
    cat("There is no result for iteration ",i,"\n")
  }
}

perf_validate_df <- perf_validate_list[lapply(perf_validate_list, length) > 0] %>% do.call("rbind",.)
perf_test_df <- perf_test_list[lapply(perf_test_list, length) > 0] %>% do.call("rbind",.)
perf_df <- rbind(perf_validate_df,perf_test_df)

write.table(perf_df,paste0("perf_fromClinical_",phenotype,".txt"), row.names = F, quote = F, sep = "\t")