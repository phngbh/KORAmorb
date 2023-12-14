library(tibble)
library(tidyverse)
library(xtable)
library(cowplot)
library(stringr)
library(scales)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(devtools)
library(limma)
library(corpcor)
library(caret)
source("UnMetImp.R")

preprocessing_transcriptomics <- function(
  # Preprocess transcriptomic data
  # Output: A matrix of preprocessed data 
  data = NULL, # matrix of raw data
  samples = NULL, # character vector of sample IDs to use 
  feature_annotation = NULL, #dataframe of feature annotation
  technical_df = NULL, # dataframe of technical variables to correct for
  technical_correction = TRUE
){
  
  # Filter for relevant samples
  data <- data[,samples]
  
  # Filter out bad quality probes
  keep = feature_annotation$ID[grep("good", feature_annotation$`QC comment`)]
  sex_tr = feature_annotation$ID[feature_annotation$ILMN_CHR %in% c("X","Y")]
  data.fil = data[ !rownames(data) %in% sex_tr & rownames(data) %in% keep, intersect(rownames(technical_df), colnames(data))]
  technical_df <- technical_df[colnames(data.fil),]
  
  if (technical_correction){
    if (all(c("p_amplification", "RIN", "sample_storage_time") %in% colnames(technical_df))){
      #Remove batch effect
      data_nobatch <- removeBatchEffect(x=data.fil, 
                                       batch = technical_df$p_amplification,
                                       covariates = technical_df[,c("RIN","sample_storage_time")])
      return(data_nobatch)
    } else {
      stop("Do not have technical variables to correct for")
    }
  } else {
    return(data)
  }
}

preprocessing_proteomics <- function(
    # Preprocess transcriptomic data
  # Output: A matrix of preprocessed data 
  data = NULL, # matrix of raw data
  samples = NULL, # character vector of sample IDs to use 
  sample_info = NULL, # dataframe of sample info
  protein_info = NULL # dataframe of protein info
){
  
  #Filter for relevant samples
  data <- data[samples, ]
  
  # Drop proteins with known bad quality
  drop = c("2795-23_3", "3590-8_3", "5071-3_3", "5073-30_2", "5118-74_2")
  protein_info = dplyr::filter(protein_info, ColCheck == "PASS" & !ID %in% drop)
  protein_info$Target[protein_info$Target == "14.03.2003"] = "14-3-3"
  
  # Drop samples with bad quality
  sample_info = dplyr::filter(sample_info, RowCheck == "PASS")
  
  # Extract remaining samples and proteins
  data = data[intersect(as.character(sample_info$SampID),rownames(data)), 
              intersect(protein_info$ID,colnames(data))] %>% 
    t() %>% log2()
  
  return(data)
}

preprocessing_metabolomics <- function(
    # Preprocess transcriptomic data
  # Output: A matrix of preprocessed data 
  data = NULL, # matrix of raw data
  samples = NULL, # character vector of sample IDs to use 
  imputation = c("load_imputed", "do"), # imputation manner 
  imputed_file = NULL # path to imputed file (optional, provide only when imputation == "load_imputed")
){
  
  #Filter for relevant samples
  data <- data[samples, ]
  
  # Filter out metabolites with > 70% NAs
  keep <- apply(data, 2, function(x) sum(!is.na(x)) >= nrow(data)/100*30)
  data <- data[,keep]
  
  if (imputation == "do"){ # Do imputation with KNN
    data.woc <- mutate(data, outcome = rep(1, times = nrow(data)))
    met <- colnames(data)
    data.imp.knn <- UnMetImp(DataFrame = data.woc, imp_type = 'knn', 
                           group1 = met, outcome = "outcome", use_covars = FALSE , logScale = F)
    data.imp <- data.imp.knn$mids %>% dplyr::select(-outcome)
    saveRDS(data.imp,"data_metab.imp.rds")
  } else {
    data.imp <- readRDS(imputed_file)
  }
  
  # Filter out outliers
  mah.dist = mahalanobis(data.imp,
                         center = colMeans(data.imp),
                         cov = pseudoinverse(cov(data.imp)),
                         inverted = T)
  mah.dist.std = abs((mah.dist - mean(mah.dist))/sd(mah.dist))
  keep = names(mah.dist.std[mah.dist.std<=4])
  data.imp = data.imp[keep,] %>% t()
  
  return(data.imp)
}

preprocessing_metabolomicsBioc <- function(
  # Preprocess metabolomic data
  # Output: A matrix of preprocessed data 
  data = NULL, # matrix of raw data
  samples = NULL, # character vector of sample IDs to use 
  technical_df = NULL # dataframe of technical variables to correct for
){
  #Extract samples with complete information
  sampls <- intersect(rownames(technical_df),intersect(rownames(data), samples))
  
  #Filter for relevant samples
  data <- data[sampls, ]
  
  #Remove technical effect
  if (!is.null(technical_df)){
    technical_df <- technical_df[sampls,,drop = F]
    data_nobatch <- removeBatchEffect(x=t(data), 
                                      batch = technical_df[,1])
    return(data_nobatch)
  } else {
    return(t(data))
  }
}

preprocessing_clinical <- function(
  # Preprocess clinical data
  # Output: A matrix of preprocessed data 
  data = NULL, # dataframe of raw data
  samples = NULL, # character vector of sample IDs to use 
  discrete_vars = NULL, # character vector of discrete variable names
  continuous_vars = NULL # character vector of continuous variable names
){
  
  #Filter for relevant samples
  data <- data[as.character(data$zz_nr) %in% samples, ]
  cat("Dimensions of the table after filtering for relevant samples: ", dim(data), "\n")
  
  # Extract discrete and continuous variables
  discrete_vars <- intersect(discrete_vars, colnames(data))
  for (v in discrete_vars){
    data[,v] <- as.character(data[,v])
    data[,v] <- ifelse(is.na(data[,v]),"missing",data[,v])
  }
  n_levels <- data %>% dplyr::select(discrete_vars) %>% apply(.,2, function(x) nlevels(as.factor(x)))
  discrete_vars <- discrete_vars[n_levels > 1]
  cat <- data %>% dplyr::select(discrete_vars) %>%
    lapply(., as.factor) %>%
    as.data.frame() 
  # cat_dum <- caret::dummyVars(~ ., data = cat)
  # cat <- as.data.frame(predict(cat_dum, newdata = cat))
  f <- reformulate(termlabels = colnames(cat), intercept = F)
  cat <- model.matrix(f, data = cat)
  continuous_vars <- intersect(continuous_vars, colnames(data))
  con <- data %>% dplyr::select(continuous_vars)
  selected_data <- cbind(cat, con)
  rownames(selected_data) <- data$zz_nr
  
  # Filter out highly missing vars
  keep <- apply(selected_data, 2, function(x) sum(!is.na(x)) > nrow(selected_data)/10*9.8)
  selected_data = selected_data[,keep]  %>% na.omit()
  cat("Dimensions of the table after filtering out highly missing variables: ", dim(selected_data), "\n")
  
  return(selected_data)
}

preprocessing_methylomics <- function(
    # Preprocess methylomic data
  # Output: A matrix of preprocessed data 
  data = NULL, # matrix of beta values
  samples = NULL, # character vector of sample IDs to use 
  technical_vars = NULL, # dataframe of technical variables
  detection_p = NULL # # dataframe of detection p-values 
){
  
  #CpG site and sample filtering
  beta.fil <- champ.filter(beta = data, pd = technical_vars, detP = detection_p, autoimpute = T, ProbeCutoff = 0.05)
  
  #Technical effect correction
  beta.fil.df = beta.fil$beta
  covariates = technical_vars[,c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")]
  beta.nobatch = removeBatchEffect(x = log(beta.fil.df),
                                   batch = info[,"Plate"],
                                   batch2 = info[,"Chip"],
                                   covariates = covariates) %>% exp()
  champ.SVD(beta = beta.nobatch, pd = beta.fil$pd) #Double-check the technical effects after correction
  
  return(beta.nobatch)
}

preprocessing_OLINK <- function(
  # Preprocess OLINK data
  # Output: A matrix of preprocessed data 
  data = NULL, # dataframe of raw data
  samples = NULL # character vector of sample IDs to use
){
  
  #Filter for relevant samples
  data <- data[rownames(data) %in% samples, ]
  
  # Filter out highly missing vars
  keep <- apply(data, 2, function(x) sum(!is.na(x)) >= nrow(data)/10*3)
  data = data[,keep]  %>% na.omit() 
  
  return(data)
}