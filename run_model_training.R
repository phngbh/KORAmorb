#!/usr/bin/env Rscript

# libs <- c("Rlib", .libPaths())
# .libPaths(libs)

cat("Load functions and environment\n")
source("model_training.R")

cat("Get the argument settings\n")
opt <- get_args()
wdir <- opt$wdir
iter <- opt$iter
outcome_name <- opt$outcome
integration <- opt$integration
algorithm <- opt$algorithm
p_metric <- opt$p_metric
features <- opt$feature
outdir <- opt$outdir

setwd(wdir)

tmp_dict <- list(Methylomics = "meth", Transcriptomics = "tra", Proteomics = "pro", MetabolomicsMetb = "metaMetb", MetabolomicsBioc = "metaBioc", Olink = "oli", Clinical = "cli")

cat("Load preprocessed data and necessary files\n")
targets <- readRDS("targets.rds")
modals_ids <- readRDS("modals_ids.rds")
preprocessed_data <- list()
for (i in names(tmp_dict)){
  if (file.exists(paste0("processed_",tmp_dict[[i]],"_",outcome_name,".rds"))){
    preprocessed_data[[i]] <- readRDS(paste0("processed_",tmp_dict[[i]],"_",outcome_name,".rds"))
    if (!i %in% c("Olink", "Clinical")){
      preprocessed_data[[i]] <- t(preprocessed_data[[i]])
    }
  }
}
# processed_meth <- readRDS(paste0("processed_meth_",outcome_name,".rds")) 
# processed_tra <- readRDS(paste0("processed_tra_",outcome_name,".rds"))
# processed_pro <- readRDS(paste0("processed_pro_",outcome_name,".rds"))
# processed_metaMetb <- readRDS(paste0("processed_metaMetb_",outcome_name,".rds"))
# processed_metaBioc <- readRDS(paste0("processed_metaBioc_",outcome_name,".rds"))
# processed_oli <- readRDS(paste0("processed_oli_",outcome_name,".rds"))
# processed_clin <- readRDS(paste0("processed_clin_",outcome_name,".rds"))
# preprocessed_data <- list(Methylomics = t(processed_meth), Transcriptomics = t(processed_tra), Proteomics = t(processed_pro), MetabolomicsMetb = t(processed_metaMetb), MetabolomicsBioc = t(processed_metaBioc), Olink = processed_oli, Clinical = processed_clin)
cat("...There are ",length(preprocessed_data), " modality loaded\n")

cat("Load feature selection results\n")
selection_result <- list()
for (i in names(tmp_dict)){
  if (file.exists(paste0("selectionRes_",tmp_dict[[i]],"_",outcome_name,".rds"))){
    selection_result[[i]] <- readRDS(paste0("selectionRes_",tmp_dict[[i]],"_",outcome_name,".rds"))
  }
}
# selectionRes_meth <- readRDS("selectionRes_meth.rds")
# selectionRes_tra <- readRDS("selectionRes_tra.rds")
# selectionRes_pro <- readRDS("selectionRes_pro.rds")
# selectionRes_metaMetb <- readRDS("selectionRes_metaMetb.rds")
# selectionRes_metaBioc <- readRDS("selectionRes_metaBioc.rds")
# selectionRes_oli <- readRDS("selectionRes_oli.rds")
# selectionRes_cli <- readRDS("selectionRes_cli.rds")
# selection_result <- list(Methylomics = selectionRes_meth, Transcriptomics = selectionRes_tra, Proteomics = selectionRes_pro, MetabolomicsMetb = selectionRes_metaMetb, MetabolomicsBioc = selectionRes_metaBioc, Olink = selectionRes_oli, Clinical = selectionRes_cli)
cat("...There are ",length(selection_result), " modality loaded\n")
if (any(!names(preprocessed_data) %in% names(selection_result)) | any(!names(selection_result) %in% names(preprocessed_data))){
  stop("The modality in the provided data and feature selection results do not match")
} else {
  preprocessed_data <- preprocessed_data[names(selection_result)]
}

cat("Make selected data list\n")
selected_data <- make_selected_list(data_list = preprocessed_data, feature_selection_result = selection_result, feature_type = features)
cat("...Free up some memory\n")
rm(preprocessed_data)
gc()
cat("...Add genomics data to selected data list\n")
if (file.exists(paste0(outcome_name,"_genomics_selected.bed"))){
  selected_data$Genomics <- bed_to_df(paste0(outcome_name,"_genomics_selected.bed"))
} else {
  cat("......There is no selected genomics data\n")
}

cat("Make input data list\n")
outcome <- targets[,outcome_name]
names(outcome) <- targets$Samples
outcome <- na.omit(outcome)
cat("There are ", length(outcome), " samples with outcome information in total\n")
input_data <- make_input_list(data_list = selected_data, outcome = outcome, id_table = modals_ids)
saveRDS(input_data,paste0(outcome_name,"_input_data.rds"))

cat("Make cross validation list\n")
outcome <- outcome[rownames(input_data$Clinical)]
cv_list <- make_cv_list(outcome = as.factor(outcome))

cat("Train and evaluate models\n")
outcome <- ifelse(outcome == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
if (integration == "FFS"){
  res <- fit_forwardSelect(data_list = input_data, y = outcome, cv_list = cv_list, p_metric = p_metric, algorithm = algorithm, n = iter)
  saveRDS(res,file.path(outdir, paste0(paste(outcome_name,integration,algorithm,p_metric,features,iter,sep = "_"),".rds")))
  res <- fit_forwardSelectFromClinical(data_list = input_data, y = outcome, cv_list = cv_list, p_metric = p_metric, algorithm = algorithm, n = iter)
  saveRDS(res,file.path(outdir, paste0(paste(outcome_name,integration,algorithm,p_metric,features,iter,"fromClinical",sep = "_"),".rds")))
} else if (integration == "ensemble") {
  res <- fit_ensemble(data_list = input_data, y = outcome, cv_list = cv_list, p_metric = p_metric, algorithm = algorithm, n = iter)
  saveRDS(res,file.path(outdir, paste0(paste(outcome_name,integration,algorithm,p_metric,features,iter,sep = "_"),".rds")))
} else {
  res <- fit_forwardSelect(data_list = input_data, y = outcome, cv_list = cv_list, p_metric = p_metric, algorithm = algorithm, n = iter)
  saveRDS(res,file.path(outdir, paste0(paste(outcome_name,"FFS",algorithm,p_metric,features,iter,sep = "_"),".rds")))
  res <- fit_forwardSelectFromClinical(data_list = input_data, y = outcome, cv_list = cv_list, p_metric = p_metric, algorithm = algorithm, n = iter)
  saveRDS(res,file.path(outdir, paste0(paste(outcome_name,"FFS",algorithm,p_metric,features,iter,"fromClinical",sep = "_"),".rds")))
  res <- fit_ensemble(data_list = input_data, y = outcome, cv_list = cv_list, p_metric = p_metric, algorithm = algorithm, n = iter)
  saveRDS(res,file.path(outdir, paste0(paste(outcome_name,"ensemble",algorithm,p_metric,features,iter,sep = "_"),".rds")))
}
