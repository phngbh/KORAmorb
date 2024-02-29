#remotes::install_github("cran/DMwR")
library(tidyverse)
library(snpStats)
library(caret)
library(caretEnsemble)
library(haven)
library(limma)
library(glmnet)
library(pROC)
library(glmnetUtils)
library(DMwR)
library(MLmetrics)
library(randomForest)
library(xgboost)
library(mlrMBO)
library(kernlab)
library(nnls)
library(caTools)
library(optparse)

get_args <- function(
  # Parse arguments passed to the R script
){
  # Define the command line options
  option_list <- list(
    make_option(c("--wdir"), type="character", default=NULL, help="Working directory"),
    make_option(c("--iter"), type="numeric", default=NULL, help="Resample/iteration"),
    make_option(c("--outcome"), type="character", default=NULL, help="Name of the outcome variable"),
    make_option(c("--integration"), type="character", default="FFS", help="Type of modality integration during model training ('FFS', 'ensemble' or 'both')"),
    make_option(c("--algorithm"), type="character", default="glmnet", help="Machine learning algorithm"),
    make_option(c("--p_metric"), type="character", default="AUPRC", help="Performance metrics"),
    make_option(c("--feature"), type="character", default="union", help="Type of features used for training ('GSEA', 'thresholding' or 'union')"),
    make_option(c("--outdir"), type="character", default=NULL, help="Output directory")
  )
  
  # Parse the command line options
  opt_parser <- OptionParser(option_list=option_list)
  opt <- parse_args(opt_parser)
  return(opt)
}

standardise <- function(
  ###Scale datasets into standard distribution
  df = NULL #dataframe or matrix to be transformed
){
  rm <- apply(as.data.frame(df), 2, function(x) sd(x) == 0)
  df <- df[,!rm, drop = FALSE]
  df <- apply(as.data.frame(df), 2, function(x) (x - mean(x))/sd(x)) %>% as.matrix()
  return(df)
}

bed_to_df = function(
  ###Load in plink bed file and make df
  bed = NULL #path to bed file
){
  
  dat <- read.plink(bed = bed)
  df <- as.data.frame(dat$genotypes)
  rn <- rownames(df)
  df.char <- apply(df,2,as.numeric)
  rownames(df.char) <- rn
  for (i in 1:nrow(df.char)){
    for (j in 1:ncol(df.char)){
      if (df.char[i,j] == 1){
        df.char[i,j] = 2
      } else if (df.char[i,j] == 2){
        df.char[i,j] = 1
      } else if (df.char[i,j] == 3){
        df.char[i,j] = 0
      }
    }
  }
  return(df.char)
}

make_selected_list <- function(
  ### Make a list of tables, each from a modality with selected features
  data_list = NULL, # list of preprocessed modality tables of shape (n_samples, n_features)
  feature_selection_result = NULL, # list of selected features
  feature_type = NULL # type of selected features
){
  selected_data <- list()
  for (m in names(feature_selection_result)){
    fs_res <- feature_selection_result[[m]]
    
    if (m == "Clinical"){
      selected_features <- fs_res$selected_vars %>% na.omit()
      
    } else {
      
      if (m %in% c("Methylomics", "Transcriptomics","Proteomics")){
        
        if (fs_res$resampling == T){
          gsea_probe <- fs_res$gsea_probe_rra %>% unlist()
          de_probe <- fs_res$d_probe_rra %>% unlist()
        } else {
          gsea_probe <- fs_res$gsea_probe %>% unlist()
          de_probe <- fs_res$d_probe %>% unlist()
        }
        
        if (feature_type == "union"){
          selected_features <- union(gsea_probe, de_probe) %>% na.omit()
        } else if (feature_type == "GSEA"){
          selected_features <- gsea_probe %>% na.omit()
        } else if (feature_type == "thresholding"){
          selected_features <- de_probe %>% na.omit()
        }
      } 
      
      if (m %in% c("MetabolomicsBioc","MetabolomicsMetb", "Olink")){
        
        if (fs_res$resampling == T){
          gsea_probe <- fs_res$msea_m_rra
          de_probe <- fs_res$d_m_rra
        } else {
          gsea_probe <- fs_res$msea_m
          de_probe <- fs_res$d_m
        }
        
        if (feature_type == "union"){
          selected_features <- union(gsea_probe, de_probe) %>% na.omit()
        } else if (feature_type == "GSEA"){
          selected_features <- gsea_probe %>% na.omit()
        } else if (feature_type == "thresholding"){
          selected_features <- de_probe %>% na.omit()
        }
      }
    }
   
    if (length(selected_features) == 0 | is.null(selected_features)){
      cat("...There is no selected features for ",m,"\n")
      selected_data[[m]] <- NULL
    } else {
      selected_data[[m]] <- data_list[[m]][,intersect(selected_features,colnames(data_list[[m]])), drop = F] 
      cat("...Dimensions of selected ", m, "table is ", dim(selected_data[[m]]),"\n")
    }
  }
  return(selected_data)
}

make_input_list <- function(
  ### Make a list modality input for training (tables of (n_samples, n_features) with matching rownames) 
  data_list = NULL, # list of modality tables of selected data
  outcome = NULL, # named character vector of outcome
  id_table = NULL # table of matching IDs, of shape (n_samples, n_modality) 
){
  # Make sure all samples have outcome info
  # id_table <- dplyr::select(id_table, all_of(names(data_list)))
  id_table <- apply(id_table,2,function(x) as.character(x)) %>% as.data.frame()
  id_table <- dplyr::filter(id_table, Clinical %in% names(outcome)) %>% na.omit()
  cat("The ID table of samples with available outcome information has dimensions ", dim(id_table),"\n")
  
  # Make sure all samples present in data list
  cat("Further filter for samples that are present in the provided data tables\n")
  id_list <- as.list(id_table)
  for (m in names(data_list)){
    id_m <- ifelse(id_list[[m]] %in% rownames(data_list[[m]]), id_list[[m]], NA)
    cat("...",m, " has ", length(na.omit(id_m)), " samples that have outcome information and are present in its corresponding data table\n")
    id_list[[m]] <- id_m
  }
  
  # Make final input list
  id_table <- as.data.frame(id_list) %>% na.omit()
  cat("The final ID table has dimensions ", dim(id_table),"\n")
  id_list <- as.list(id_table)
  for (m in names(data_list)){
    data_list[[m]] <- data_list[[m]][id_list[[m]],,drop = F] %>% standardise()
    rownames(data_list[[m]]) <- id_table$Clinical
  }
  return(data_list)
}

make_cv_list <- function(
  ###Create a list of cross-validation subsets, same for all feature groups
  outcome = NULL, #named vector of outcome (factor) 
  times = 100, # number of repeated resamples
  partition = 0.8, # fraction of train set per resample
  kfold_inner = 5, #how many inner folds
  times_inner = 1, #number of k-fold cross validations for inner training
  seed = 993 #setting fixed seeds 
){
  
  set.seed(seed)
  outer_train <- createDataPartition(y = outcome, times = times, p = partition)
  outer_test <- lapply(outer_train, function(x) setdiff(1:length(outcome), x))
  cat("The outer train sets have length ", length(outer_train[[1]]), "\n")
  cat("The outer test sets have length ", length(outer_test[[1]]), "\n")
  
  outer <- list(train = lapply(outer_train, function(x) names(outcome)[x]), test = lapply(outer_test, function(x) names(outcome)[x]))
  
  inner_train <- lapply(outer$train, function(x) {set.seed(seed); createMultiFolds(y = outcome[x], k = kfold_inner, times = times_inner)})
  cat("The inner train sets have length ", length(inner_train[[1]][[1]]), "\n")
  outcome_outer <- lapply(outer$train, function(x) names(outcome[x]))
  inner_test <- list()
  for (i in 1:length(inner_train)){
    inner_test[[i]] <- list()
    for ( j in 1: length(inner_train[[i]])){
      inner_train[[i]][[j]] <- outcome_outer[[i]][inner_train[[i]][[j]]]
      inner_test[[i]][[j]] <- setdiff(outcome_outer[[i]], inner_train[[i]][[j]])
    }
  }
  cat("The inner test sets have length ", length(inner_test[[1]][[1]]), "\n")
  inner <- list(train = inner_train, test = inner_test)
  return(list(outer = outer, inner = inner))
}

LogLoss <- function(pred, true, eps = 1e-15, weights = NULL) {

  pred = pmin(pmax(pred, eps), 1 - eps) # Bound the results
  
  if (is.null(weights)) {
    return(-(sum(
      true * log(pred) + (1 - true) * log(1 - pred)
    )) / length(true))
  } else{
    return(-weighted.mean(true * log(pred) + (1 - true) * log(1 - pred), weights))
  }
}

caretLogLoss <- function(data, lev = NULL, model = NULL) {
  cls <- levels(data$obs) #find class names
  loss <- LogLoss(
    pred = data[, cls[2]],
    true = as.numeric(data$obs) - 1,
    weights = data$weights
  )
  names(loss) <- c('myLogLoss')
  loss
}

bestSE = function(
  #Customize selection function for parametter tunning (mean + 1SE) 
  x, metric, num, maximize
){
  
  perf <- x[, metric] + (x[, paste(metric, "SD", sep = "")])/sqrt(num)
  if (maximize){
    bestIter = which.max(perf)
  } else {
    bestIter = which.min(perf)
  }
  bestIter
}

null_model <- function(
  # Make a random prediction considering class imbalance
  target = NULL, # character vector/ factor of target variable
  positive = NULL # character string of positive (minority) class
){
  n <- length(target)
  negative <- setdiff(levels(as.factor(target)), positive)
  # Determine the proportion of positive class examples
  positive_proportion <- sum(target == positive) / n
  
  # Generate random predictions with the same proportion
  random_predictions <- sample(c(negative, positive), n, replace = TRUE, prob = c(1 - positive_proportion, positive_proportion))
  
  return(random_predictions)
}

train_model <- function(
  ### Train a model and return prediction performance, prediction probabilities (and optionally important variables)
  x_train = NULL, # matrix of shape (n_samples,n_features)
  y_train = NULL, # factor of two levels in order Case & Control
  x_test = NULL, # matrix of shape (n_samples,n_features)
  y_test = NULL, # factor of two levels in order Case & Control
  algorithm = c("glmnet","rf","svmRadial","xgbTree"),
  p_metric = c("wLogLoss","AUROC","AUPRC"),
  return_imp_vars = TRUE,
  seed = 993
){
  
  if (p_metric == "wLogLoss"){
    sampling <- NULL
    sumFunc <- caretLogLoss
    metric <- "myLogLoss"
    maximize <- F
    weights_train <- ifelse(y_train == "One", table(y_train)[[2]]/table(y_train)[[1]], 1)
    weights_test <- ifelse(y_test == "One", table(y_test)[[2]]/table(y_test)[[1]], 1)
  } else if (p_metric == "AUROC"){
    sampling <- "smote"
    sumFunc <- twoClassSummary
    metric <- "ROC"
    maximize <- T
    weights_train <- NULL
    weights_test <- ifelse(y_test == "One", table(y_test)[[2]]/table(y_test)[[1]], 1)
  } else if (p_metric == "AUPRC"){
    sampling <- "smote"
    sumFunc <- prSummary
    metric <- "AUC"
    maximize <- T
    weights_train <- NULL
    weights_test <- ifelse(y_test == "One", table(y_test)[[2]]/table(y_test)[[1]], 1)
  }
  
  my_control <- trainControl(
    method="repeatedcv",
    number=5,
    repeats = 4,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=sumFunc,
    sampling = sampling,
    allowParallel = T,
    verboseIter = F
  )
  
  # if (algorithm == "xgbTree"){
  #   grid <- expand.grid(
  #     nrounds = c(10, 50, 100),           # Number of boosting rounds
  #     max_depth = c(3, 6, 9),              # Maximum depth of a tree
  #     eta = c(0.01, 0.1, 0.3),             # Learning rate
  #     gamma = c(0, 0.2, 0.4),              # Minimum loss reduction for a split
  #     colsample_bytree = c(0.6, 0.8, 1.0),  # Subsample ratio of columns
  #     min_child_weight = seq(1, 50,5),
  #     subsample = seq(0.1, 1,0.1)
  #   )
  # } else {
  #   grid <- NULL
  # }
  
  set.seed(seed)
  fit <- caret::train(x = x_train,
                      y = y_train,
                      method= algorithm, 
                      metric=metric,
                      #tuneGrid = grid,
                      tuneLength = 20,
                      weights = weights_train, 
                      maximize = maximize,
                      trControl=my_control,
                      importance = TRUE,
                      verbosity = 0)
  pred = predict(fit, x_test, s = "lambda.min", type = "prob")$One
  names(pred) <- names(y_test)
  roc <- roc(response = y_test, predictor = pred, levels = c("Zero","One"))
  roc = auc(roc)[[1]]
  pr = MLmetrics::PRAUC(pred, ifelse(y_test == "One",1,0))
  ll = LogLoss(pred, ifelse(y_test == "One",1,0), weights = weights_test)
  
  if (return_imp_vars){
    var_imp <- varImp(fit)$importance
    if (algorithm == "glmnet"){
      var <- var_imp$Overall
    } else {
      var <- var_imp$One
    }
    names(var) <- rownames(var_imp)
    var <- var[order(var, decreasing = T)]
  } else {
    var <- NULL
  }
  
  return(list(roc = roc, pr = pr, ll = ll, pred = pred, var = var))
  
}

fit_forwardSelect <- function(
  #Function to do CV and identify best model through forward feature selection
  data_list = NULL, # List of input table each per modality (matched rownames)
  y = NULL, #Named vector of outcome
  cv_list = NULL, # List generated by make_cv_list()
  p_metric = c("wLogLoss","AUROC","AUPRC"),
  algorithm = c("glmnet","rf","svmRadial"), #Core learning algorithm for the FFS
  seed = 993,
  n = NULL #Iteration number in parallelizing 
){
  
  i <- n
  cat("Start process for iteration ",i,"\n")
  
  y_train_i <- y[cv_list$outer$train[[i]]]
  cat("There are ",length(y_train_i), " samples in train set\n")
  y_test_i <- y[cv_list$outer$test[[i]]]
  cat("There are ",length(y_test_i), " samples in test set\n")
  n_cv <- length(cv_list$inner$train[[i]]) 
  cat("Number of inner CV folds is ",n_cv,"\n")
  
  if (nlevels(as.factor(y_train_i)) < 2 | nlevels(as.factor(y_test_i)) < 2){
    stop("Train or test set has only one class label!")
  }
  
  cat("Set up training parameters\n")
  cat("...Algorithm is ",algorithm,"\n")
  cat("...Optimization metric is ",p_metric,"\n")
  
  cat("Make all combinations of modality\n")
  n_datasets <- length(data_list)
  cat("There are ",n_datasets," modality to analyze\n")
  comb_list <- lapply(1:n_datasets, function(x) combn(names(data_list),x, simplify = F))
  
  perf_validate <- list()
  perf_test <- list()
  pred_list <- list()
  var_list <- list()
  best_d <- NULL
  
  cat("Iterate through combinations\n")
  for (j in 1:length(comb_list)){
    
    cat("Evaluate combinations of ",j, " modality\n")
    perf_j <- list(roc = list(), pr = list(), ll = list())
    
    #Filter for combinations that involve best performing modality in the previously smaller combination 
    if (!is.null(best_d)){
      ind <- lapply(comb_list[[j]], function(x) all(best_d %in% x)) %>% unlist()
      comb_list_fil <- comb_list[[j]][ind]
    } else {
      comb_list_fil <- comb_list[[j]]
    }
    
    for (d in 1:length(comb_list_fil)){
      
      comb <- paste0(comb_list_fil[[d]], collapse = "")
      cat("...Evaluate ",comb,"\n")
      
      roc_d <- list()
      pr_d <- list()
      ll_d <- list()
      for (n in 1:n_cv){
        cat("......Fold ",n,"\n")
        x_train_n <- do.call(cbind, data_list[comb_list_fil[[d]]])[cv_list$inner$train[[i]][[n]],, drop = F] %>% as.data.frame()
        if (ncol(x_train_n) == 1){
          x_train_n <- cbind(x_train_n, ranv = 0)
        }
        y_train_n <- y[rownames(x_train_n) ]
        x_test_n <- do.call(cbind, data_list[comb_list_fil[[d]]])[cv_list$inner$test[[i]][[n]],, drop = F] %>% as.data.frame()
        if (ncol(x_test_n) == 1){
          x_test_n <- cbind(x_test_n, ranv = 0)
        }
        y_test_n <- y[rownames(x_test_n) ]
        
        res <- train_model(x_train = x_train_n, y_train = y_train_n, x_test = x_test_n, y_test = y_test_n, algorithm = algorithm, p_metric = p_metric, return_imp_vars = F, seed = seed)
        roc_d[[n]] = res$roc 
        pr_d[[n]] = res$pr 
        ll_d[[n]] = res$ll 
      }
      
      perf_j$roc[[d]] <- unlist(roc_d) %>% median()
      perf_j$pr[[d]] <- unlist(pr_d) %>% median()
      perf_j$ll[[d]] <- unlist(ll_d) %>% median()
      cat("......Median prediction performance on test sets: AUROC: ",perf_j$roc[[d]],", AUPRC: ", perf_j$pr[[d]], " & weighted log loss: ",perf_j$ll[[d]],"\n")
      
    }
    
    if (p_metric == "AUROC"){
      best_ind <- which.max(unlist(perf_j$roc))
    } else if (p_metric == "AUPRC") {
      best_ind <- which.max(unlist(perf_j$pr))
    } else if (p_metric == "wLogLoss") {
      best_ind <- which.min(unlist(perf_j$ll))
    }
    
    best_d <- comb_list_fil[[best_ind]]
    cat("...The best ",j, "-modality model is ",best_d,"\n")
    perf_validate[[j]] <- data.frame(Complexity = paste0(j,"_modality"), Model = paste0(best_d, collapse = ""), Value = c(perf_j$roc[[best_ind]],perf_j$pr[[best_ind]],perf_j$ll[[best_ind]]), Type = c("AUROC","AUCPR","Weighted LogLoss"))
    
    cat("...Evaluating on outer test set\n")
    x_train_i <- do.call(cbind, data_list[comb_list_fil[[best_ind]]])[cv_list$outer$train[[i]],, drop = F] %>% as.data.frame()
    if (ncol(x_train_i) == 1){
      x_train_i <- cbind(x_train_i, ranv = 0)
    }
    x_test_i <- do.call(cbind, data_list[comb_list_fil[[best_ind]]])[cv_list$outer$test[[i]],, drop = F] %>% as.data.frame()
    if (ncol(x_test_i) == 1){
      x_test_i <- cbind(x_test_i, ranv = 0)
    }
    
    if (any(rownames(x_train_i) != names(y_train_i)) | any(rownames(x_test_i) != names(y_test_i))){
      stop("Samples in train and test sets do not match!")
    }
    
    res_test <- train_model(x_train = x_train_i, y_train = y_train_i, x_test = x_test_i, y_test = y_test_i, algorithm = algorithm, p_metric = p_metric, return_imp_vars = T, seed = seed)
    pred_test <- res_test$pred
    roc_test <- res_test$roc
    pr_test <- res_test$pr
    ll_test <- res_test$ll
    var <- res_test$var
    
    pred_list[[j]] <- pred_test
    var_list[[j]] <- var
    
    cat("......Prediction performance on outer test set: AUROC: ",roc_test,", AUPRC: ", pr_test, " & weighted log loss: ",ll_test,"\n")
    
    perf_test[[j]] <- data.frame(Complexity = paste0(j,"_modality"), Model = paste0(best_d, collapse = ""), Value = c(roc_test,pr_test,ll_test), Type = c("AUROC","AUCPR","Weighted LogLoss"))
    
  }
  
  perf_validate <- do.call(rbind, perf_validate)
  perf_test <- do.call(rbind, perf_test)
  cat("Done\n")
  return(list(perf_validate = perf_validate, perf_test = perf_test, var = var_list, predProbs_test = pred_list, test_samples = cv_list$outer$test[[i]]))
}

fit_forwardSelectFromClinical = function(
  #Function to do CV and identify best model through foward feature selection
  data_list = NULL,
  y = NULL, #Named vector of outcome
  cv_list = NULL,
  p_metric = c("wLogLoss","AUROC","AUPRC"),
  algorithm = c("glmnet","rf","svmRadial"),
  seed = 993,
  n = NULL #Fold number in parallelizing
){
  
  i <- n
  
  cat("Start process for iteration ",i,"\n")
  
  y_train_i <- y[cv_list$outer$train[[i]]]
  cat("There are ",length(y_train_i), " samples in train set\n")
  y_test_i <- y[cv_list$outer$test[[i]]]
  cat("There are ",length(y_test_i), " samples in test set\n")
  n_cv <- length(cv_list$inner$train[[i]]) 
  cat("Number of inner CV folds is ",n_cv,"\n")
  
  if (nlevels(as.factor(y_train_i)) < 2 | nlevels(as.factor(y_test_i)) < 2){
    stop("Train or test set has only one class label!")
  }
  
  cat("Set up training parameters\n")
  cat("...Algorithm is ",algorithm,"\n")
  cat("...Optimization metric is ",p_metric,"\n")
  
  cat("Cross-validation with clinical as baseline\n")
  roc_clin = list()
  pr_clin = list()
  ll_clin = list()
  for (n in 1:n_cv){
    cat("...Fold ",n,"\n")
    x_train_n_clin = data_list$Clinical[cv_list$inner$train[[i]][[n]],] %>% as.data.frame()
    if (ncol(x_train_n_clin) == 1){
      x_train_n_clin <- cbind(x_train_n_clin, ranv = 0)
    }
    y_train_n_clin = y[rownames(x_train_n_clin)]
    x_test_n_clin = data_list$Clinical[cv_list$inner$test[[i]][[n]],] %>% as.data.frame()
    if (ncol(x_test_n_clin) == 1){
      x_test_n_clin <- cbind(x_test_n_clin, ranv = 0)
    }
    y_test_n_clin = y[rownames(x_test_n_clin)]
    
    res <- train_model(x_train = x_train_n_clin, y_train = y_train_n_clin, x_test = x_test_n_clin, y_test = y_test_n_clin, algorithm = algorithm, p_metric = p_metric, return_imp_vars = F, seed = seed)
    roc_clin[[n]] = res$roc 
    pr_clin[[n]] = res$pr 
    ll_clin[[n]] = res$ll 
  }
  
  cat("...Median prediction performance on test sets: AUROC: ",median(unlist(roc_clin)),", AUPRC: ", median(unlist(pr_clin)), " & weighted log loss: ",median(unlist(ll_clin)),"\n")
  perf_clin_validate = data.frame(Complexity = "1_modality", Model = "Clinical", Value = c(median(unlist(roc_clin)), median(unlist(pr_clin)), median(unlist(ll_clin))), Type = c("AUROC","AUPRC","Weighted LogLoss"))
  
  cat("...Testing with clinical as baseline\n")
  x_train_i_clin = data_list$Clinical[cv_list$outer$train[[i]],] %>% as.data.frame()
  if (ncol(x_train_i_clin) == 1){
    x_train_i_clin <- cbind(x_train_i_clin, ranv = 0)
  }
  x_test_i_clin = data_list$Clinical[cv_list$outer$test[[i]],] %>% as.data.frame()
  if (ncol(x_test_i_clin) == 1){
    x_test_i_clin <- cbind(x_test_i_clin, ranv = 0)
  }
  
  if (any(rownames(x_train_i_clin) != names(y_train_i)) | any(rownames(x_test_i_clin) != names(y_test_i))){
    stop("Samples in train and test sets do not match!")
  }
  
  res_test <- train_model(x_train = x_train_i_clin, y_train = y_train_i, x_test = x_test_i_clin, y_test = y_test_i, algorithm = algorithm, p_metric = p_metric, return_imp_vars = T, seed = seed)
  pred_test_clin <- res_test$pred
  roc_test_clin <- res_test$roc
  pr_test_clin <- res_test$pr
  ll_test_clin <- res_test$ll
  var_clin <- res_test$var
  
  cat("...Prediction performance on outer test set: AUROC: ",roc_test_clin,", AUPRC: ", pr_test_clin, " & weighted log loss: ",ll_test_clin,"\n")
  perf_clin_test = data.frame(Complexity = "1_modality", Model = "Clinical", Value = c(roc_test_clin,pr_test_clin,ll_test_clin), Type = c("AUROC","AUPRC","Weighted LogLoss"))
  
  cat("Make all combinations out of the rest of modality\n")
  n_datasets <- length(data_list)
  cat("There are ",n_datasets - 1," modality to analyze\n")
  comb_list = lapply(1:(n_datasets - 1), function(x) combn(setdiff(names(data_list),"Clinical"),x, simplify = F))
  
  perf_validate = list()
  perf_test = list()
  pred_list <- list()
  var_list = list()
  best_d = NULL
  
  cat("Iterate through combinations\n")
  for (j in 1:length(comb_list)){
    
    cat("Evaluate combinations of ",j, " modality\n")
    perf_j = list(roc = list(), pr = list(), ll = list())
    
    #Filter for combinations that involve best performing modality in the previously smaller combination 
    if (!is.null(best_d)){
      ind = lapply(comb_list[[j]], function(x) all(best_d %in% x)) %>% unlist()
      comb_list_fil = comb_list[[j]][ind]
    } else {
      comb_list_fil = comb_list[[j]]
    }
    
    for (d in 1:length(comb_list_fil)){
      
      comb = paste0(c("Clinical",comb_list_fil[[d]]), collapse = "")
      cat("...Evaluate ",comb,"\n")
      
      roc_d = list()
      pr_d = list()
      ll_d = list()
      for (n in 1:n_cv){
        
        cat("......Fold ",n,"\n")
        x_train_n = do.call(cbind, data_list[c("Clinical",comb_list_fil[[d]])])[cv_list$inner$train[[i]][[n]],] %>% as.data.frame()
        if (ncol(x_train_n) == 1){
          x_train_n <- cbind(x_train_n, ranv = 0)
        }
        y_train_n = y[rownames(x_train_n) ]
        x_test_n = do.call(cbind, data_list[c("Clinical",comb_list_fil[[d]])])[cv_list$inner$test[[i]][[n]],] %>% as.data.frame()
        if (ncol(x_test_n) == 1){
          x_test_n <- cbind(x_test_n, ranv = 0)
        }
        y_test_n = y[rownames(x_test_n) ]
        
        res <- train_model(x_train = x_train_n, y_train = y_train_n, x_test = x_test_n, y_test = y_test_n, algorithm = algorithm, p_metric = p_metric, return_imp_vars = F, seed = seed)
        roc_d[[n]] = res$roc 
        pr_d[[n]] = res$pr 
        ll_d[[n]] = res$ll 
      }
      perf_j$roc[[d]] = unlist(roc_d) %>% median()
      perf_j$pr[[d]] = unlist(pr_d) %>% median()
      perf_j$ll[[d]] = unlist(ll_d) %>% median()
    }
    
    if (p_metric == "AUROC"){
      best_ind = which.max(unlist(perf_j$roc))
    } else if (p_metric == "AUPRC") {
      best_ind <- which.max(unlist(perf_j$pr))
    } else if (p_metric == "wLogLoss") {
      best_ind = which.min(unlist(perf_j$ll))
    }
    
    best_d = comb_list_fil[[best_ind]]
    cat("...The best ",j+1, "-modality model is ",paste0(c("Clinical",best_d), collapse = ""),"\n")
    
    perf_validate[[j]] = data.frame(Complexity = paste0(j+1,"_modality"), Model = paste0(c("Clinical",best_d), collapse = ""), Value = c(perf_j$roc[[best_ind]],perf_j$pr[[best_ind]],perf_j$ll[[best_ind]]), Type = c("AUROC","AUPRC","Weighted LogLoss"))
    
    cat("...Evaluating on outer test set\n")
    x_train_i = do.call(cbind, data_list[comb_list_fil[[best_ind]]])[cv_list$outer$train[[i]],] %>% as.data.frame()
    if (ncol(x_train_i) == 1){
      x_train_i <- cbind(x_train_i, ranv = 0)
    }
    x_test_i = do.call(cbind, data_list[comb_list_fil[[best_ind]]])[cv_list$outer$test[[i]],] %>% as.data.frame()
    if (ncol(x_test_i) == 1){
      x_test_i <- cbind(x_test_i, ranv = 0)
    }
    
    if (any(rownames(x_train_i) != names(y_train_i)) | any(rownames(x_test_i) != names(y_test_i))){
      stop("Samples in train and test sets do not match!")
    }
    
    res_test <- train_model(x_train = x_train_i, y_train = y_train_i, x_test = x_test_i, y_test = y_test_i, algorithm = algorithm, p_metric = p_metric, return_imp_vars = T, seed = seed)
    pred_test <- res_test$pred
    roc_test <- res_test$roc
    pr_test <- res_test$pr
    ll_test <- res_test$ll
    var <- res_test$var
    
    pred_list[[j]] <- pred_test
    var_list[[j]] = var
    perf_test[[j]] = data.frame(Complexity = paste0(j+1,"_modality"), Model = paste0(c("Clinical",best_d), collapse = ""), Value = c(roc_test,pr_test,ll_test), Type = c("AUROC","AUPRC","Weighted LogLoss"))
    cat("......Prediction performance on outer test set: AUROC: ",roc_test,", AUPRC: ", pr_test, " & weighted log loss: ",ll_test,"\n")
  }
  
  perf_validate = do.call(rbind, perf_validate) %>% rbind(perf_clin_validate,.)
  perf_test = do.call(rbind, perf_test) %>% rbind(perf_clin_test,.)
  cat("Done\n")
  return(list(perf_validate = perf_validate, perf_test = perf_test, var = c(var_clin,var_list), predProbs_test = c(pred_test_clin, pred_list), train_samples = cv_list$outer$train[[i]])) 
}

fit_ensemble = function(
  data_list = NULL,
  y = NULL, #Named vector of outcome
  algorithm = c("glmnet","rf","svmRadial","xbgTree"),
  p_metric = c("AUROC","AUPRC","wLogLoss"),
  cv_list = NULL,
  seed = 993,
  n = NULL #Fold number in parallelizing 
){
  
  df = do.call(cbind, data_list) %>% as.data.frame()
  ef = list()
  i = n
    
  cat("Iteration ",i,"\n")
  
  x_train_i = df[cv_list$outer$train[[i]],,drop=F]
  if (ncol(x_train_i) == 1){
    x_train_i <- cbind(x_train_i, ranv = 0)
  }
  x_test_i = df[cv_list$outer$test[[i]],,drop=F]
  if (ncol(x_test_i) == 1){
    x_test_i <- cbind(x_test_i, ranv = 0)
  }
  
  y_train = y[cv_list$outer$train[[i]]]
  y_test = y[cv_list$outer$test[[i]]]
  
  weights = ifelse(y_train == "One", table(y_train)[[2]]/table(y_train)[[1]], 1)
  weights_test = ifelse(y_test == "One", table(y_test)[[2]]/table(y_test)[[1]], 1)
  
  if (nlevels(as.factor(y_train)) < 2 | nlevels(as.factor(y_test)) < 2){
    stop("Train or test set has only one class label!")
  }
  
  if (any(rownames(x_train_i) != names(y_train)) | any(rownames(x_test_i) != names(y_test))){
    stop("Samples in train and test sets do not match!")
  }
  
  if (p_metric == "wLogLoss"){
    sampling <- NULL
    sumFunc <- caretLogLoss
    metric <- "myLogLoss"
    maximize <- F
    weights_train <- ifelse(y_train == "One", table(y_train)[[2]]/table(y_train)[[1]], 1)
    weights_test <- ifelse(y_test == "One", table(y_test)[[2]]/table(y_test)[[1]], 1)
  } else if (p_metric == "AUROC"){
    sampling <- "smote"
    sumFunc <- twoClassSummary
    metric <- "ROC"
    maximize <- T
    weights_train <- NULL
    weights_test <- ifelse(y_test == "One", table(y_test)[[2]]/table(y_test)[[1]], 1)
  } else if (p_metric == "AUPRC"){
    sampling <- "smote"
    sumFunc <- prSummary
    metric <- "AUC"
    maximize <- T
    weights_train <- NULL
    weights_test <- ifelse(y_test == "One", table(y_test)[[2]]/table(y_test)[[1]], 1)
  }
  
  my_control <- trainControl(
    method="repeatedcv",
    number=5,
    repeats = 4,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=sumFunc,
    sampling = sampling,
    allowParallel = T,
    verboseIter = F
  )
  
  cat("Fit concatened dataset\n")
  set.seed(seed)
  fit_cat <- caret::train(x = x_train_i,
                          y = y_train,
                          method=algorithm, 
                          metric=metric,
                          tuneLength = 20,
                          weights = weights_train, 
                          maximize = maximize,
                          trControl=my_control,
                          importance = TRUE)
  
  pred_cat = predict(fit_cat, x_test_i, s = "lambda.min", type = "prob")
  
  cat("Fit individual datasets\n")
  #ef[[i]] = list()
  
  sigl_var <- F 
  for (d in names(data_list)){
    
    cat("...Data type: ",d,"\n")
    
    x_train_d = data_list[[d]][cv_list$outer$train[[i]],,drop=F]
    if (ncol(x_train_d) == 1){
      x_train_d <- cbind(x_train_d, ranv = 0)
      sigl_var <- T
    }
    #x_test = data_list[[d]][cv_list$outer$test[[i]],]
    
    set.seed(seed)
    fit_d <- caret::train(x = x_train_d,
                        y = y_train,
                        method=algorithm, 
                        metric=p_metric,
                        tuneLength = 20,
                        weights = weights_train, 
                        maximize = maximize,
                        trControl=my_control,
                        importance = TRUE)
    ef[[d]] = fit_d
    
  }
  
  cat("Fit meta model and computing performance\n")
  preds = list()
  vars = list()
  for (d in names(ef)){
    preds[[d]] = arrange(ef[[d]]$pred, Resample, rowIndex)$One
    vars[[d]] = varImp(ef[[d]])$importance #%>% mutate(Importance = (.$One + .$Zero)/2) %>% dplyr::select(-One,-Zero) 
    rownames(vars[[d]]) = rownames(varImp(ef[[d]])$importance)
  }
  preds = as.data.frame(preds)
  y_train_truth = y_train[arrange(ef[[1]]$pred, Resample, rowIndex)$rowIndex] 
  y_train_truth = ifelse(y_train_truth == "One", 1, 0)
  y_train_best = abs(preds - y_train_truth) %>% apply(.,1,which.min)
  y_train_best = colnames(preds)[y_train_best]
  cat("...Best datasets include:\n")
  print(table(y_train_best))
  set.seed(seed)
  fit_meta = caret::train(#x = preds,
                          #y = y_train[arrange(ef[[i]][[1]]$pred, Resample, rowIndex)$rowIndex],
                          x = x_train_i[arrange(ef[[1]]$pred, Resample, rowIndex)$rowIndex,,drop = F],
                          y = y_train_best %>% as.factor(),
                          method="glmnet", 
                          metric="AUC",
                          tuneLength = 20, 
                          trControl=trainControl(
                            method="boot",
                            number=10,
                            savePredictions="final",
                            classProbs=TRUE,
                            summaryFunction=multiClassSummary
                          ),
                          importance = TRUE)
  
  if (sigl_var){
    x_test_i <- cbind(x_test_i, ranv = 0)
  }
  pred_mods <- predict(ef, x_test_i, type = "prob") %>% as.data.frame()
  keep_col = grep("One", colnames(pred_mods))
  pred_mods = pred_mods[,keep_col, drop = F]
  colnames(pred_mods) = names(ef)
  # pred_ensemble = predict(fit_meta, pred_mods, s = "lambda.min", type = "prob")
  # pred_mods$ensemble = pred_ensemble$One
  pred_weight = predict(fit_meta, x_test_i, s = "lambda.min", type = "prob")
  pred_mods_fil = pred_mods[,colnames(pred_weight)]
  # pred_best = predict(fit_meta, x_test_i, s = "lambda.min", type = "raw") %>% as.character()
  # pred_mods$ensemble = lapply(1:nrow(pred_mods), function (x) pred_mods[x,pred_best[x]]) %>% unlist()
  pred_mods$concatenation = pred_cat$One
  pred_mods$ensemble = apply(pred_mods_fil*pred_weight,1,sum)
  # y_test = ifelse(y_test == "One",1,0)
  perf_roc = caTools::colAUC(pred_mods, y_test)
  perf_pr = lapply(pred_mods, function(x) MLmetrics::PRAUC(x, ifelse(y_test == "One",1,0))) %>% as.data.frame()
  perf_ll = lapply(pred_mods, function(x) LogLoss(x, ifelse(y_test == "One",1,0), weights = weights_test)) %>% as.data.frame()
  # cf = coef(fit_meta$finalModel, fit_meta$bestTune$mtry)
  # weight = cf[-1]
  # names(weight) = names(vars)
  # vars_df = lapply(names(vars), function(x) vars[[x]]*abs(weight[x])) %>% do.call(rbind, .)
  # var_imp = varImp(fit_meta)$importance #%>% mutate(Importance = (.$One + .$Zero)/2) %>% dplyr::select(-One,-Zero)
  # rownames(var_imp) = names(vars)
  # vars_df = lapply(names(vars), function(x) vars[[x]]*abs(var_imp[x,1])) %>% do.call(rbind, .)
  # vars_df$Importance = vars_df$Importance/sum(vars_df$Importance)*100
  # vars_df$Overall = vars_df$Overall/sum(vars_df$Overall)*100
  
  # cat("...Meta model and computing performance\n")
  # class(ef[[i]]) <- c("caretList")
  # x_test_i = df[cv_list$outer$test[[i]],]
  # pred_mods <- predict(ef[[i]], x_test_i) %>% as.data.frame()
  # colnames(pred_mods) = names(ef[[i]])
  # model_ensemble <- caretStack(
  #   ef[[i]],
  #   method = "glmnet",
  #   metric="ROC",
  #   trControl=trainControl(
  #     method = "boot",
  #     number=10,
  #     savePredictions="final",
  #     summaryFunction=twoClassSummary,
  #     classProbs=TRUE
  #   ))
  # pred_ensemble = predict(model_ensemble$ens_model, pred_mods, type = "prob")
  # pred_mods$ensemble = pred_ensemble$One 
  # perf = caTools::colAUC(pred_mods, y_test)
  # var_imp = varImp(model_ensemble)
  
  res <- list(mod_list = ef ,
                 mod_meta = fit_meta,
                 pred_mods = pred_mods,
                 perf = list(AUROC = perf_roc, AUPRC = perf_pr, wLogLoss = perf_ll),
                 vars = vars,
                 weight = pred_weight
                 #coef = cf
                 )
  
  cat("Done\n")
  return(res)
}

fit_forwardSelectEnsemble = function(
  #Function to do CV and identify best model through foward feature selection
  data_list = NULL,
  y = NULL, #Named vector of outcome
  cv_list = NULL,
  n_datasets = NULL, #Number of datasets to evaluate
  metric = c("wLogLoss","AUROC"),
  seed = 993,
  n = NULL #Fold number in parallelizing 
){
  
  df = do.call(cbind, data_list) %>% as.data.frame()
  i = n
  
  cat("Folds ",i,"\n")
  
  x_train_i = df[cv_list$outer$train[[i]],]
  x_test_i = df[cv_list$outer$test[[i]],]
  
  y_train = y[cv_list$outer$train[[i]]]
  y_test = y[cv_list$outer$test[[i]]]
  
  if (nlevels(as.factor(y_train)) < 2 | nlevels(as.factor(y_test)) < 2){
    stop("Train or test set has only one class label!")
  }
  
  if (any(rownames(x_train_i) != names(y_train)) | any(rownames(x_test_i) != names(y_test))){
    stop("Samples in train and test sets do not match!")
  }
  
  if (metric == "wLogLoss"){
    sampling = NULL
    sumFunc = caretLogLoss
    metric = "myLogLoss"
    maximize = F
    weights = ifelse(y_train == "One", table(y_train)[[2]]/table(y_train)[[1]], 1)
    weights_test = ifelse(y_test == "One", table(y_test)[[2]]/table(y_test)[[1]], 1)
  } else if (metric == "AUROC"){
    sampling = "smote"
    sumFunc = twoClassSummary
    metric = "ROC"
    maximize = T
    weights = rep(1, length(y_train))
    weights_test = ifelse(y_test == "One", table(y_test)[[2]]/table(y_test)[[1]], 1)
  }
  
  my_control <- trainControl(
    method="repeatedcv",
    #number=5,
    #repeats = 4,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=sumFunc,
    sampling = sampling,
    index = lapply(cv_list$inner_train[[i]],function(x) match(x, names(y_train))),
    allowParallel = T
  )
  
  cat("Fit base models\n")
  base = list()
  for (d in names(data_list)){
    
    cat("...Data type: ",d,"\n")
    
    x_train_d = data_list[[d]][cv_list$outer$train[[i]],]
    
    set.seed(seed)
    fit_d <- caret::train(x = x_train_d,
                          y = y_train,
                          method="glmnet", 
                          metric=metric,
                          tuneLength = 20,
                          weights = weights, 
                          maximize = maximize,
                          trControl=my_control,
                          importance = TRUE)
    base[[d]] = fit_d
  }
  
  cat("Evaluating base models performances\n")
  pred_base_test <- predict(base, x_test_i, type = "prob") %>% as.data.frame()
  keep_col = grep("One", colnames(pred_base_test))
  pred_base_test = pred_base_test[,keep_col, drop = F]
  colnames(pred_base_test) = names(base)
  roc_base = caTools::colAUC(pred_base_test, ifelse(y_test == "One",1,0))
  ll_base = lapply(pred_base_test, function(x) LogLoss(x, ifelse(y_test == "One",1,0), weights = weights_test)) %>% as.data.frame()
  pr_base = lapply(pred_base_test, function(x) MLmetrics::PRAUC(x, ifelse(y_test == "One",1,0))) %>% as.data.frame()
  best_ind_roc_base = apply(roc_base,1,which.max)
  best_ind_ll_base = apply(ll_base,1,which.min)
  best_ind_pr_base = apply(pr_base,1,which.max)
  best_base = colnames(pred_base_test)[best_ind_roc_base]
  cat("...The best base model is ",best_base,"\n")
  
  return(list(perf = data.frame(Complexity = "1 Dataset", 
                                Model = c(colnames(roc_base)[best_ind_roc_base],colnames(pr_base)[best_ind_pr_base], colnames(ll_base)[best_ind_ll_base]), 
                                Value = c(roc_base[1,best_ind_roc_base],pr_base[1,best_ind_pr_base],ll_base[1,best_ind_ll_base]), 
                                Type = c("AUROC","AUCPR","Weighted LogLoss")), mod = base[[best_base]]))
  
  cat("Prepare for ensemble learning\n")
  pred_base_train = lapply(base, function(x) arrange(x$pred, Resample, rowIndex)$One) %>% as.data.frame()
  y_train_truth = y_train[arrange(base[[1]]$pred, Resample, rowIndex)$rowIndex]
  y_train_truth = ifelse(y_train_truth == "One", 1, 0)
  comb_list = lapply(2:n_datasets, function(x) combn(names(data_list),x, simplify = F))

  perf_all = list()
  fit_all = list()
  best_c = NULL

  for (j in 1:length(comb_list)){

    ind_j = lapply(comb_list[[j]], function(x) best_base %in% x) %>% unlist()
    comb_list_j = comb_list[[j]][ind_j]

    perf = list(roc = list(), pr = list(), ll = list())
    fit = list()

    if (!is.null(best_c)){
      ind = lapply(comb_list_j, function(x) all(best_c %in% x)) %>% unlist()
      comb_list_fil = comb_list_j[ind]
    } else {
      comb_list_fil = comb_list_j
    }

    for (c in 1:length(comb_list_fil)){

      comb = paste0(comb_list_fil[[c]], collapse = "")
      cat("Fit ensemble of ",comb,"\n")

      pred_base_c = pred_base_train[,comb_list_fil[[c]]]
      y_train_best_ind = abs(pred_base_c - y_train_truth) %>% apply(.,1,which.min)
      y_train_best = colnames(pred_base_c)[y_train_best_ind]
      # cat("...Best datasets include:\n")
      # print(table(y_train_best))

      set.seed(seed)
      fit_meta = caret::train(x = x_train_i[arrange(base[[1]]$pred, Resample, rowIndex)$rowIndex,],
                              y = y_train_best %>% as.factor(),
                              method="glmnet",
                              metric="logLoss",
                              tuneLength = 20,
                              trControl=trainControl(
                                method="boot",
                                number=10,
                                savePredictions="final",
                                classProbs=TRUE,
                                summaryFunction=multiClassSummary
                              ),
                              importance = TRUE)

      pred_weight = predict(fit_meta, x_test_i, s = "lambda.min", type = "prob")
      pred_base_fil = pred_base_test[,colnames(pred_weight)]
      pred_ensemble = apply(pred_base_fil*pred_weight,1,sum)
      perf_roc = caTools::colAUC(pred_ensemble, y_test)[,1]
      perf_pr = MLmetrics::PRAUC(pred_ensemble, ifelse(y_test == "One",1,0))
      perf_ll = LogLoss(pred_ensemble, ifelse(y_test == "One",1,0), weights = weights_test)

      fit[[comb]] = fit_meta
      perf$roc[[comb]] = perf_roc
      perf$pr[[comb]] = perf_pr
      perf$ll[[comb]] = perf_ll
    }

    best_ind_roc = unlist(perf$roc) %>% which.max()
    best_ind_pr = unlist(perf$pr) %>% which.max()
    best_ind_ll = unlist(perf$ll) %>% which.min()
    perf_all[[j]] = data.frame(Complexity = paste0(j+1," Datasets"), Model = c(names(fit)[best_ind_roc],names(fit)[best_ind_pr], names(fit)[best_ind_ll]),
                               Value = c(perf$roc[[best_ind_roc]], perf$pr[[best_ind_pr]], perf$ll[[best_ind_ll]]), Type = c("AUROC","AUCPR","Weighted LogLoss"))
    fit_all[[j]] = fit[[best_ind_ll]]

    best_c = comb_list_fil[[best_ind_ll]]
    cat("The best ensemble of",best_base,"and", j, "base model is",best_c,"\n")
  }

  perf_all = do.call(rbind, perf_all) %>% rbind(data.frame(Complexity = "1 Dataset", Model = c(colnames(roc_base)[best_ind_roc_base],colnames(pr_base)[best_ind_pr_base], colnames(ll_base)[best_ind_ll_base]),
                                                           Value = c(roc_base[1,best_ind_roc_base],pr_base[1,best_ind_pr_base],ll_base[1,best_ind_ll_base]), Type = c("AUROC","AUCPR","Weighted LogLoss")),.)
  return(list(perf = perf_all, mod = fit_all))
  
}

benchmark <- function(
  # Benchmark performance of our final model with null and/or published models
  df = NULL,
  target = NULL,
  positive = NULL, # character string of positive /minority class
  final_vars = NULL,
  published_vars = NULL,
  p_metric = c("wLogLoss","AUROC","AUPRC"),
  algorithm = c("glmnet","rf","svmRadial"), #Core learning algorithm for the FFS
  n_iterations = 100,
  seed = 993
){
  cat("Initiate result list and make resamples\n")
  res <- list(auc_fin = list(), auc_ran = list(), auc_est = list())
  set.seed(seed)
  resamples <- createDataPartition(y = target, times = n_iterations, p = 0.8)
  
  cat("Iterate through resamples\n")
  for (i in 1:length(resamples)){
    cat("Iter ",i,"\n")
    train_samples <- rownames(df)[resamples[[i]]]
    x_train <- df[train_samples, , drop = FALSE]
    y_train <- target[train_samples]
    x_test <- df[setdiff(rownames(df), train_samples),,drop = FALSE]
    y_test <- target[rownames(x_test)]
    
    cat("...Make random prediction\n")
    pred_ran <- null_model(target = y_test, positive = positive) 
    pred_ran <- ifelse(pred_ran == positive,1,0)
    res$auc_ran[[i]] = MLmetrics::PRAUC(pred_ran, ifelse(y_test == positive,1,0))
    cat("......AUPRC: ",res$auc_ran[[i]],"\n")
    
    # cat("...Make final model prediction\n")
    # x_train_fin <- x_train[,final_vars,drop=FALSE]
    # if (ncol(x_train_fin) == 1){
    #   x_train_fin <- cbind(x_train_fin, ranv = 0)
    # }
    # x_test_fin <- x_test[,final_vars,drop=FALSE]
    # if (ncol(x_test_fin) == 1){
    #   x_test_fin <- cbind(x_test_fin, ranv = 0)
    # }
    # res_fin <- train_model(x_train = x_train_fin, y_train = y_train, x_test = x_test_fin, y_test = y_test, algorithm = algorithm, p_metric = p_metric, return_imp_vars = F, seed = seed)
    # res$auc_fin[[i]] <- res_fin$pr
    # cat("......AUPRC: ",res$auc_fin[[i]],"\n")
    
    if(!is.null(published_vars)){
      cat("...Make published model prediction\n")
      x_train_est <- x_train[,published_vars,drop=FALSE]
      if (ncol(x_train_est) == 1){
        x_train_est <- cbind(x_train_est, ranv = 0)
      }
      x_test_est <- x_test[,published_vars,drop=FALSE]
      if (ncol(x_test_est) == 1){
        x_test_est <- cbind(x_test_est, ranv = 0)
      }
      res_est <- train_model(x_train = x_train_est, y_train = y_train, x_test = x_test_est, y_test = y_test, algorithm = algorithm, p_metric = p_metric, return_imp_vars = F, seed = seed)
      res$auc_est[[i]] <- res_est$pr
      cat("......AUPRC: ",res$auc_est[[i]],"\n")
    }
  }
  
  cat("Done\n")
  return(res)
  
}