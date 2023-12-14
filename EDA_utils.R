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
library(caret)

make_stratified_samplings <- function(
  # Create resamples of a set of samples that retain the class distribution of the target
  # Output: A list of resamples
  samples = NULL, # a character vector of sample IDs
  target = NULL, # a named character vector of binary target (0=ctrl, 1=case)
  p = 0.8, # fraction of samples to select for each resample
  times = 100 # number of resamples
){
  # Extract relevant samples
  int <- intersect(samples, names(target))
  target <- target[int] #%>% as.character()
  
  # Create list of resamples
  sample_list <- caret::createDataPartition(y = target, times = times, p = p, list = TRUE)
  final_sample_list <- lapply(sample_list, function(x) names(target[x])) 
  
  return(final_sample_list)
}

make_CKD_target <- function(
  # Create CKD target from raw variables using pre-defined formular
  # Output: A sample ID-named binary vector for incident CKD (0 = case, 1 = control) 
  data = NULL, # dataframe containing raw variables
  id_col = "", #string specifying the column name of sample IDs
  target_vars = c("egfr_time0" = "", "egfr_time1" = "", "uacr_time0" = "", "uacr_time1" = "") # Named string vector of variable names required for calculation
  
){
  # Extract required variables
  egfr_time0 <- target_vars["egfr_time0"]
  egfr_time1 <- target_vars["egfr_time1"]
  uacr_time0 <- target_vars["uacr_time0"]
  uacr_time1 <- target_vars["uacr_time1"]
  
  # Determine control at baseline
  baseline_ctrl <- data[id_col][data[egfr_time0] >= 60 & data[uacr_time0] < 30] %>% na.omit() %>% as.character() # eGFR >= 60 ml/min/1.73m2 & uACR < 30 mg/g
  
  # Determine incident control
  ctrl <- data[id_col][data[egfr_time1] >= 60 & data[uacr_time1] < 30] %>% na.omit() %>% as.character()
  ctrl <- intersect(ctrl, baseline_ctrl)
  
  # Determine incident case
  case <- data[id_col][data[egfr_time1] < 60 | data[uacr_time1] >= 30] %>% na.omit() %>% as.character()
  case <- intersect(case, baseline_ctrl)
  
  target <- c(rep(0, length(ctrl)), rep(1, length(case)))
  names(target) <- c(ctrl, case)
  
  return(target)
}

make_DSPN_target <- function(
  # Create DSPN target from raw variables using pre-defined formular
  # Output: A sample ID-named binary vector for incident DSPN (0 = case, 1 = control) 
  data = NULL, # dataframe containing raw variables
  id_col = "", #string specifying the column name of sample IDs
  target_vars = c("mnsi_time0" = "", "mnsi_time1" = "") # Named string vector of variable names required for calculation
  
){
  # Extract required variables
  mnsi_time0 <- target_vars["mnsi_time0"]
  mnsi_time1 <- target_vars["mnsi_time1"]
  
  # Determine control at baseline
  baseline_ctrl <- data[id_col][data[mnsi_time0] < 3] %>% na.omit() %>% as.character() # MNSI score < 3
  
  # Determine incident control
  ctrl <- data[id_col][data[mnsi_time1] < 3] %>% na.omit() %>% as.character()
  ctrl <- intersect(ctrl, baseline_ctrl)
  
  # Determine incident case
  case <- data[id_col][data[mnsi_time1] >= 3] %>% na.omit() %>% as.character()
  case <- intersect(case, baseline_ctrl)
  
  target <- c(rep(0, length(ctrl)), rep(1, length(case)))
  names(target) <- c(ctrl, case)
  
  return(target)
}

make_sleepDisorder_target <- function(
  # Create sleep disorder target from raw variables using pre-defined formular
  # Output: A sample ID-named binary vector for incident sleep disorder (0 = case, 1 = control) 
  data = NULL, # dataframe containing raw variables
  id_col = "", #string specifying the column name of sample IDs
  target_vars = c("sleep1_time0" = "", "sleep1_time1" = "", "sleep2_time0" = "", "sleep2_time1" = "") # Named string vector of variable names required for calculation
  
){
  # Extract required variables
  sleep1_time0 <- target_vars["sleep1_time0"]
  sleep1_time1 <- target_vars["sleep1_time1"]
  sleep2_time0 <- target_vars["sleep2_time0"]
  sleep2_time1 <- target_vars["sleep2_time1"]
  
  # Determine control at baseline
  baseline_ctrl <- data[id_col][data[sleep1_time0] == 3 & data[sleep2_time0] == 3] %>% na.omit() %>% as.character() # Never had sleep problem
  
  # Determine incident control
  ctrl <- data[id_col][data[sleep1_time1] == 3 & data[sleep2_time1] == 3] %>% na.omit() %>% as.character()
  ctrl <- intersect(ctrl, baseline_ctrl)
  
  # Determine incident case
  case <- data[id_col][data[sleep1_time1] != 3 | data[sleep2_time1] != 3] %>% na.omit() %>% as.character()
  case <- intersect(case, baseline_ctrl)
  
  target <- c(rep(0, length(ctrl)), rep(1, length(case)))
  names(target) <- c(ctrl, case)
  
  return(target)
}

make_hypertension_target <- function(
  # Create hypertension target from raw variables using pre-defined formular
  # Output: A sample ID-named binary vector for incident hypertension (0 = case, 1 = control) 
  data = NULL, # dataframe containing raw variables
  id_col = "", #string specifying the column name of sample IDs
  target_vars = c("hyper_time0" = "", "hyper_time1" = "") # Named string vector of variable names required for calculation
  
){
  # Extract required variables
  hyper_time0 <- target_vars["hyper_time0"]
  hyper_time1 <- target_vars["hyper_time1"]
  
  # Determine control at baseline
  baseline_ctrl <- data[id_col][data[hyper_time0] == 2 ] %>% na.omit() %>% as.character() 
  
  # Determine incident control
  ctrl <- data[id_col][data[hyper_time1] == 2] %>% na.omit() %>% as.character()
  ctrl <- intersect(ctrl, baseline_ctrl)
  
  # Determine incident case
  case <- data[id_col][data[hyper_time1] == 1] %>% na.omit() %>% as.character()
  case <- intersect(case, baseline_ctrl)
  
  target <- c(rep(0, length(ctrl)), rep(1, length(case)))
  names(target) <- c(ctrl, case)
  
  return(target)
}

make_dyslipidemia_target <- function(
  # Create dyslipidemia target from raw variables using pre-defined formular
  # Output: A sample ID-named binary vector for incident dyslipidemia (0 = case, 1 = control) 
  data = NULL, # dataframe containing raw variables
  id_col = "", #string specifying the column name of sample IDs
  target_vars = c("ldl_time0" = "", "ldl_time1" = "", 
                  #"hdl_time0" = "", "hdl_time1" = "",
                  #"chol_time0" = "", "chol_time1" = "", 
                  "trig_time0" = "", "trig_time1" = "") # Named string vector of variable names required for calculation
  
){
  # Extract required variables
  ldl_time0 <- target_vars["ldl_time0"]
  ldl_time1 <- target_vars["ldl_time1"]
  # hdl_time0 <- target_vars["hdl_time0"]
  # hdl_time1 <- target_vars["hdl_time1"]
  # chol_time0 <- target_vars["chol_time0"]
  # chol_time1 <- target_vars["chol_time1"]
  trig_time0 <- target_vars["trig_time0"]
  trig_time1 <- target_vars["trig_time1"]
  
  # Determine control at baseline
  baseline_ctrl <- data[id_col][data[ldl_time0] < 160 & #data[hdl_time0] >= 40 & data[chol_time0] < 240 & 
                                  data[trig_time0] <= 200] %>% 
    na.omit() %>% as.character() 
  
  # Determine incident control
  ctrl <- data[id_col][data[ldl_time1] < 160 & #data[hdl_time1] >= 40 & data[chol_time1] < 240 & 
                         data[trig_time1] <= 200] %>% 
    na.omit() %>% as.character()
  ctrl <- intersect(ctrl, baseline_ctrl)
  
  # Determine incident case
  case <- data[id_col][data[ldl_time1] >= 160 | #data[hdl_time1] < 40 | data[chol_time1] >= 240 | 
                         data[trig_time1] <= 200] %>% 
    na.omit() %>% as.character()
  case <- intersect(case, baseline_ctrl)
  
  target <- c(rep(0, length(ctrl)), rep(1, length(case)))
  names(target) <- c(ctrl, case)
  
  return(target)
}

make_diabetes_target <- function(
  # Create diabetes target from raw variables using pre-defined formular
  # Output: A sample ID-named binary vector for incident diabetes (0 = case, 1 = control) 
  data = NULL, # dataframe containing raw variables
  id_col = "", #string specifying the column name of sample IDs
  target_vars = c("dm_time0" = "", "dm_time1" = "") # Named string vector of variable names required for calculation
  
){
  # Extract required variables
  dm_time0 <- target_vars["dm_time0"]
  dm_time1 <- target_vars["dm_time1"]
  
  # Determine control at baseline
  baseline_ctrl <- data[id_col][data[dm_time0] == 0 ] %>% na.omit() %>% as.character() # Normal patients
  
  # Determine incident control
  ctrl <- data[id_col][data[dm_time1] == 0] %>% na.omit() %>% as.character()
  ctrl <- intersect(ctrl, baseline_ctrl)
  
  # Determine incident case
  case <- data[id_col][data[dm_time1] != 0] %>% na.omit() %>% as.character()
  case <- intersect(case, baseline_ctrl)
  
  target <- c(rep(0, length(ctrl)), rep(1, length(case)))
  names(target) <- c(ctrl, case)
  
  return(target)
}

make_MI_target <- function(
  # Create myocardial infarction target from raw variables using pre-defined formular
  # Output: A sample ID-named binary vector for incident MI (0 = case, 1 = control) 
  data = NULL, # dataframe containing raw variables
  id_col = "", #string specifying the column name of sample IDs
  target_vars = c("MI_time0" = "", "MI_time1" = "") # Named string vector of variable names required for calculation
  
){
  # Extract required variables
  MI_time0 <- target_vars["MI_time0"]
  MI_time1 <- target_vars["MI_time1"]
  
  # Determine control at baseline
  baseline_ctrl <- data[id_col][data[MI_time0] == 0] %>% na.omit() %>% as.character() 
  
  # Determine incident control
  ctrl <- data[id_col][data[MI_time1] == 0] %>% na.omit() %>% as.character()
  ctrl <- intersect(ctrl, baseline_ctrl)
  
  # Determine incident case
  case <- data[id_col][data[MI_time1] == 1] %>% na.omit() %>% as.character()
  case <- intersect(case, baseline_ctrl)
  
  target <- c(rep(0, length(ctrl)), rep(1, length(case)))
  names(target) <- c(ctrl, case)
  
  return(target)
}

make_stroke_target <- function(
    # Create myocardial infarction target from raw variables using pre-defined formular
  # Output: A sample ID-named binary vector for incident MI (0 = case, 1 = control) 
  data = NULL, # dataframe containing raw variables
  id_col = "", #string specifying the column name of sample IDs
  target_vars = c("stroke_time0" = "", "stroke_time1" = "") # Named string vector of variable names required for calculation
  
){
  # Extract required variables
  stroke_time0 <- target_vars["stroke_time0"]
  stroke_time1 <- target_vars["stroke_time1"]
  
  # Determine control at baseline
  baseline_ctrl <- data[id_col][data[stroke_time0] == 0] %>% na.omit() %>% as.character() 
  
  # Determine incident control
  ctrl <- data[id_col][data[stroke_time1] == 0] %>% na.omit() %>% as.character()
  ctrl <- intersect(ctrl, baseline_ctrl)
  
  # Determine incident case
  case <- data[id_col][data[stroke_time1] == 1] %>% na.omit() %>% as.character()
  case <- intersect(case, baseline_ctrl)
  
  target <- c(rep(0, length(ctrl)), rep(1, length(case)))
  names(target) <- c(ctrl, case)
  
  return(target)
}

make_targets <- function(
  #Create targets from raw variables based on pre-defined formular 
  #Output: A data frame of (n_samples, n_targets) and plots (optional)
  data = NULL, # dataframe containing raw variables
  id_col = "", #string specifying the column name of sample IDs
  targets = c("CKD", "DSPN", "sleep disorder", "hypertension", "dyslipidemia", "diabetes", "MI", "stroke"), # string vector of target names
  target_var_list = list("CKD" = c("egfr_time0" = "utgfr_ckd_crcc", "egfr_time1" = "u3tgfr_ckd_crcc", "uacr_time0" = "utacr", "uacr_time1" = "u3tacr"),
                          "DSPN" = c("mnsi_time0" = "utmnsi", "mnsi_time1" = "u3tmnsi"),
                          "sleep disorder" = c("sleep1_time0" = "uc099", "sleep1_time1" = "u3c099", "sleep2_time0" = "uc100", "sleep2_time1" = "u3c100"),
                          "dyslipidemia" = c("ldl_time0" = "ul_ldla", "ldl_time1" = "u3lk_ldla", 
                                             # "hdl_time0" = "ul_hdla", "hdl_time1" = "u3lk_hdla",
                                             # "chol_time0" = "ul_chola", "chol_time1" = "u3lk_chola", 
                                             "trig_time0" = "ul_tria", "trig_time1" = "u3lk_tria"),
                          "hypertension" = c("hyper_time0" = "uthyact", "hyper_time1" = "u3thyact"),
                          "diabetes" = c("dm_time0" = "utdm100y15", "dm_time1" = "u3tdm100y15"),
                         "MI" = c("MI_time0" = "prev_mi_ff4", "MI_time1" = "inz_mi_ff4"),
                         "stroke" = c("stroke_time0" = "prev_apo_ff4", "stroke_time1" = "inz_apo_ff4"))
){
  
  # Make a list of target labels
  target_list <- list()
  for (i in targets){
    if (i == "CKD"){
      target_list[[i]] <- make_CKD_target(data = data, id_col = id_col, target_vars = target_var_list[[i]])
    } else if (i == "DSPN"){
      target_list[[i]] <- make_DSPN_target(data = data, id_col = id_col, target_vars = target_var_list[[i]])
    } else if (i == "sleep disorder"){
      target_list[[i]] <- make_sleepDisorder_target(data = data, id_col = id_col, target_vars = target_var_list[[i]])
    } else if (i == "hypertension"){
      target_list[[i]] <- make_hypertension_target(data = data, id_col = id_col, target_vars = target_var_list[[i]])
    } else if (i == "dyslipidemia"){
      target_list[[i]] <- make_dyslipidemia_target(data = data, id_col = id_col, target_vars = target_var_list[[i]])
    } else if (i == "diabetes"){
      target_list[[i]] <- make_diabetes_target(data = data, id_col = id_col, target_vars = target_var_list[[i]])
    } else if (i == "MI"){
      target_list[[i]] <- make_MI_target(data = data, id_col = id_col, target_vars = target_var_list[[i]])
    } else if (i == "stroke"){
      target_list[[i]] <- make_stroke_target(data = data, id_col = id_col, target_vars = target_var_list[[i]])
    }
  }
  
  # Extract a union of samples
  samples <- lapply(target_list, function(x) names(x)) %>% unlist() %>% unique()
  
  # Rearrange the samples into same order
  target_list <- lapply(target_list, function(x) x[samples])
  # for (i in names(target_list)){
  #   ind <- match(samples, names(target_list[[i]]))
  #   target_list[[i]] <- target_list[[i]][ind]
  # }
  
  # Make final dataframe
  df <- data.frame(Samples = samples) %>% cbind(.,data.frame(target_list))
  
  return(df)
}

make_modals_list = function(
  # Make list of feature selection sets and training set for each target
  #Output: A list of sample subsets
  data = NULL, # dataframe containing raw variables
  id_col = "", #string specifying the column name of sample IDs
  #targets = c("CKD", "DSPN", "sleep disorder", "hypertension", "dyslipidemia", "diabetes"), # string vector of target names
  targets_df = NULL, # dataframe of available targets (n_samples, n_targets)
  modal_vars = c("Genomics" = "lg_dnaAxiom_s4f4", "Methylomics" = "un_meth450k_f4", "Transcriptomics" = "un_expr_f4ogtt", "Proteomics" = "un_protSoma_f4", "MetabolomicsMetb" = "un_metabMetabolon_f4", "MetabolomicsBioc" = "un_metabBiocrates_f4", "Clinical" = "zz_nr")
){
  
  # Extract modality variables
  modals_df <- dplyr::select(data, modal_vars) 
  rownames(modals_df) <- data[,id_col]
  colnames(modals_df) <- names(modal_vars)
  modals_df$Olink <- ifelse(is.na(data$uh_o_il18), NA, modals_df$Clinical)
  
  sample_list <- list()
  # Filter for target-specific samples
  for (i in colnames(targets)){
    if (i == "Samples"){
      next
    }
    sample_list[[i]] <- list()
    
    modals_df_i <- modals_df[as.character(targets_df$Samples[!is.na(targets_df[,i])]),]
    cnt <- apply(modals_df_i, 1, function(x) sum(!is.na(x))) 
    train = rownames(modals_df_i)[cnt == ncol(modals_df)] # Samples having all modalities are used for training
    sample_list[[i]][["train"]] <- train
    
    # Make feature selection sets for each modality
    fs <- modals_df_i[setdiff(rownames(modals_df_i), train),]
    sample_list[[i]][["feature_selection"]] <- list()
    for (m in colnames(modals_df_i)){
      fs_m <- na.omit(fs[,m]) %>% as.character()
      sample_list[[i]][["feature_selection"]][[m]] <- fs_m
    }
  }
  
  return(sample_list)
  
}

plot_sample_sets <- function(
  # Make barplot of subsets of samples for a list of targets
  modals_list = NULL, # List of partitioned samples
  target_df = NULL,
  modals_ids = NULL # Table of corresponding IDs for all modality
){
  
  plot_list <- list()
  for (i in names(modals_list)){
    
    plot_list[[i]] <- list()
    train_df <- data.frame(Labels = target_df[target_df$Samples %in% modals_list[[i]]$train,i]) 
    train_df$Labels <- as.factor(train_df$Labels)
    p <- ggplot(train_df, aes(x = Labels, fill = Labels)) + 
        geom_bar() + 
        theme(legend.position = "none",panel.background = element_rect(fill = "white"),
              axis.text = element_text(size = 20),
              axis.title.y = element_blank(),
              axis.title.x = element_text(size = 15)) + 
        xlab("Train") +
      coord_cartesian(ylim = c(0,nrow(train_df) + 10))
    plot_list[[i]][["train"]] <- p
    
    max_n <- lapply(modals_list[[i]]$feature_selection, function(x) length(x)) %>% unlist() %>% max()
    for (m in names(modals_list[[i]]$feature_selection)){
      
      ind <- match(modals_list[[i]]$feature_selection[[m]], as.character(modals_ids[,m]))
      samples_m_clin <- as.character(modals_ids$Clinical[ind]) %>% na.omit()
      fs_df <- data.frame(Labels = target_df[target_df$Samples %in% samples_m_clin,i]) 
      fs_df$Labels <- as.factor(fs_df$Labels)
      p <- ggplot(fs_df, aes(x = Labels, fill = Labels)) + 
        geom_bar() + 
        theme(legend.position = "none",panel.background = element_rect(fill = "white"),
              axis.text = element_text(size = 20),
              axis.title.y = element_blank(),
              axis.title.x = element_text(size = 15)) + 
        xlab(paste0("FS_",m))+
        coord_cartesian(ylim = c(0,max_n + 10))
      plot_list[[i]][[paste0("FS_",m)]] <- p
    }
  }
  
  return(plot_list)
}