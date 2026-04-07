library(bootnet)      
library(qgraph)       
library(psychonetrics) 
library(dplyr)
library(parSim)       
library(tidyr)
library(ggplot2)
library(mice)         
library(EMgaussian)   
library(tidyverse)
library(missForest)   
library(VIM)          
library(mixgb)
library(cglasso) 
library(Matrix)

# SIM

parSim(
  n_people  = c(250, 500, 750),     
  n_nodes   = 8,                    
  missing   = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),            
  estimator = c(
    "ggm_prune", 
    "EBICglasso_listwise", 
    "EBICglasso_pairwise", 
    "EBICglasso_fiml",
    "EBICglasso_single_imputation_pmm", # essentially regression-based imputation with an added matching step.
    "EBICglasso_miceMI", 
    "ggmModSelect_listwise",
    "ggmModSelect_pairwise",
    "ggmModSelect_fiml",
    "ggmModSelect_single_imputation_pmm", # essentially regression-based imputation with an added matching step.
    "ggmModSelect_miceMI",
    "EM_EBIC", 
    "EM_CV",
    # Impute inline -> EBICglasso
    "EBICglasso_mean", 
    "EBICglasso_median",
    "EBICglasso_knn", 
    "EBICglasso_missForest",
    "ggmModSelect_mean", 
    "ggmModSelect_median", 
    "ggmModSelect_knn",
    "ggmModSelect_missForest",
    "EBICglasso_xgboost"),
  
  reps   = 1,          
  write  = TRUE,       
  name   = "sims_batch_100",
  nCores = 4,         
  progressbar = FALSE, 
  
  expression = {
    
    library(bootnet); library(qgraph); library(psychonetrics)
    library(dplyr); library(mice); library(EMgaussian) 
    library(missForest); library(VIM); library(mixgb); library(cglasso)
    
    # Minimum number of imputations with a nonzero edge required for inclusion (mice)
    threshold <- 3
    
    # 1) TRUE NETWORK
    trueNetwork <- genGGM(n_nodes, propPositive = 0.8, p = 0.2)
    
    # 2) SIMULATED DATA + MISSINGNESS
    simFun  <- ggmGenerator(missing = missing)  # create random data generating function
    simData <- simFun(trueNetwork, n = n_people) # use funct to draw n obs + randomly delete
    
    # 3) ESTIMATE NETWORKS
    
    # A. Gaussian Graphical Model (Prune, EBIC + DELETION/FMIL) ---
    
    if (estimator == "ggm_prune") {
      mod <- ggm(simData, estimator = "FIML") %>% runmodel() 
      SEs <- !any(is.na(mod@parameters$se[mod@parameters$matrix=="omega"]))
      if (SEs){
        mod <- mod %>% prune(alpha = 0.05)
        estNetwork <- getmatrix(mod, "omega") 
      } else {
        stop("SEs could not be computed")
      }
      
    } else if (estimator == "EBICglasso_listwise") {
      estNetwork <- estimateNetwork(simData, default = "EBICglasso", missing = "listwise")$graph
      
    } else if (estimator == "EBICglasso_pairwise") {
      estNetwork <- estimateNetwork(simData, default = "EBICglasso", missing = "pairwise")$graph
      
    } else if (estimator == "EBICglasso_fiml") {
      modS <- varcov(simData, estimator = "FIML") %>% runmodel()
      S_fiml <- getmatrix(modS, "sigma")
      estNetwork <- qgraph::EBICglasso(S_fiml, n = nrow(simData))
      
      #  B. MICE into EBIC ---
      
    } else if (estimator == "EBICglasso_single_imputation_pmm") {
      imp <- mice(simData, m = 1, method = "pmm", maxit = 5, printFlag = FALSE)
      dat_mi <- complete(imp, 1)
      estNetwork <- estimateNetwork(dat_mi, default = "EBICglasso")$graph
      
    } else if (estimator == "EBICglasso_miceMI") {
      imp <- mice(simData, m = 5, method = "pmm", maxit = 5, printFlag = FALSE)
      graphs <- lapply(1:5, function(i) {
        dat_i <- complete(imp, i)
        estimateNetwork(dat_i, default = "EBICglasso")$graph
      })
      arr <- simplify2array(graphs)
      estNetwork <- apply(arr, c(1, 2), function(x) {
        nonzero <- x[x != 0]
        if (length(nonzero) >= threshold) median(nonzero) else 0
      }) 
      
      # C. SAME BUT FOR GGMODSELECT ---
      
    } else if (estimator == "ggmModSelect_listwise") {
      estNetwork <- estimateNetwork(simData, default = "ggmModSelect", missing = "listwise")$graph
      
    } else if (estimator == "ggmModSelect_pairwise") {
      estNetwork <- estimateNetwork(simData, default = "ggmModSelect", missing = "pairwise")$graph
      
    } else if (estimator == "ggmModSelect_fiml") {
      modS <- varcov(simData, estimator = "FIML") %>% runmodel()
      S_fiml <- getmatrix(modS, "sigma")
      estNetwork <- qgraph::ggmModSelect(S_fiml, n = nrow(simData))$graph
      
    } else if (estimator == "ggmModSelect_single_imputation_pmm") {
      imp <- mice(simData, m = 1, method = "pmm", maxit = 5, printFlag = FALSE)
      dat_mi <- complete(imp, 1)
      estNetwork <- estimateNetwork(dat_mi, default = "ggmModSelect")$graph
      
    } else if (estimator == "ggmModSelect_miceMI") {
      imp <- mice(simData, m = 5, method = "pmm", maxit = 5, printFlag = FALSE)
      graphs <- lapply(1:5, function(i) {       
        dat_i <- complete(imp, i)               
        estimateNetwork(dat_i, default = "ggmModSelect")$graph
      })
      arr <- simplify2array(graphs)          
      estNetwork <- apply(arr, c(1, 2), function(x) {
        nonzero <- x[x != 0]
        if (length(nonzero) >= threshold) median(nonzero) else 0
      })
      
      # D. Using EM with EBIC, CV ---
      
    } else if (estimator == "EM_EBIC") {
      X   <- scale(as.matrix(simData))
      rho <- EMgaussian::rhogrid(100, method = "qgraph", dat = simData)
      estNetwork <- estimateNetwork(X, fun = EMgaussian::EMggm, rho = rho,
                                    glassoversion = "glasso", rhoselect = "ebic")$graph
      
    } else if (estimator == "EM_CV") {
      X   <- scale(as.matrix(simData))
      rho <- EMgaussian::rhogrid(100, method = "qgraph", dat = X)
      estNetwork <- estimateNetwork(X, fun = EMgaussian::EMggm, rho = rho,
                                    glassoversion = "glasso", rhoselect = "kfold", k = 5)$graph
      
    } else if (estimator == "EBICglasso_mean") {
      dat_imp <- as.data.frame(simData)
      for (j in seq_along(dat_imp)) {       
        if (is.numeric(dat_imp[[j]])) {
          na_idx <- is.na(dat_imp[[j]])
          dat_imp[[j]][na_idx] <- mean(dat_imp[[j]], na.rm = TRUE)
        }
      }
      estNetwork <- estimateNetwork(dat_imp, default = "EBICglasso")$graph
      
    } else if (estimator == "EBICglasso_median") {
      dat_imp <- as.data.frame(simData)
      for (j in seq_along(dat_imp)) {
        if (is.numeric(dat_imp[[j]])) {
          na_idx <- is.na(dat_imp[[j]])
          dat_imp[[j]][na_idx] <- stats::median(dat_imp[[j]], na.rm = TRUE)
        }
      }
      estNetwork <- estimateNetwork(dat_imp, default = "EBICglasso")$graph
      
    } else if (estimator == "EBICglasso_knn") {
      dat_imp <- VIM::kNN(simData, k = 5, imp_var = FALSE)
      estNetwork <- estimateNetwork(dat_imp, default = "EBICglasso")$graph
      
    } else if (estimator == "EBICglasso_missForest") {
      dat_imp <- missForest(as.data.frame(simData), maxiter = 5, ntree = 100, verbose = FALSE)$ximp
      estNetwork <- estimateNetwork(dat_imp, default = "EBICglasso")$graph
      
    } else if (estimator == "ggmModSelect_mean") {
      dat_imp <- as.data.frame(simData)
      for (j in seq_along(dat_imp)) {
        if (is.numeric(dat_imp[[j]])) {
          na_idx <- is.na(dat_imp[[j]])
          dat_imp[[j]][na_idx] <- mean(dat_imp[[j]], na.rm = TRUE)
        }
      }
      estNetwork <- estimateNetwork(dat_imp, default = "ggmModSelect")$graph
      
    } else if (estimator == "ggmModSelect_median") {
      dat_imp <- as.data.frame(simData)
      for (j in seq_along(dat_imp)) {
        if (is.numeric(dat_imp[[j]])) {
          na_idx <- is.na(dat_imp[[j]])
          dat_imp[[j]][na_idx] <- stats::median(dat_imp[[j]], na.rm = TRUE)
        }
      }
      estNetwork <- estimateNetwork(dat_imp, default = "ggmModSelect")$graph
      
    } else if (estimator == "ggmModSelect_knn") {
      dat_imp <- VIM::kNN(simData, k = 5, imp_var = FALSE)
      estNetwork <- estimateNetwork(dat_imp, default = "ggmModSelect")$graph
      
    } else if (estimator == "ggmModSelect_missForest") {
      dat_imp <- missForest(as.data.frame(simData), maxiter = 5, ntree = 100, verbose = FALSE)$ximp
      estNetwork <- estimateNetwork(dat_imp, default = "ggmModSelect")$graph
      
      # F. XG BOOST into EBIC/GGM (using mixgb package) (Deng & Lumley) ---
      
    } else if (estimator == "EBICglasso_xgboost") {
      dat_xgb <- as.data.frame(simData)
      dat_xgb[] <- lapply(dat_xgb, as.numeric)  
      dat_imp <- mixgb(
        data = dat_xgb, m = 1, nrounds = 50, maxit = 1, verbose = FALSE,
        xgb.params = list(max_depth = 4, eta = 0.3, nthread = 1))
      dat_imp <- as.data.frame(dat_imp)
      estNetwork <- estimateNetwork(dat_imp, default = "EBICglasso")$graph
      
    } else {
      stop("Invalid estimator")
    } 
    
    
    # 4) RECOVERY METRICS + EDGE COUNTS + ERRORS
    lt <- lower.tri(trueNetwork, diag = FALSE)
    tE <- trueNetwork[lt]
    eE <- estNetwork[lt]
    
    # confusion counts
    TP <- sum(tE != 0 & eE != 0)
    FP <- sum(tE == 0 & eE != 0)
    FN <- sum(tE != 0 & eE == 0)
    TN <- sum(tE == 0 & eE == 0)
    
    # edge counts (will need this to show weird findings later)
    n_true_edges <- sum(tE != 0)
    n_estimated_edges <- sum(eE != 0)
    
    # metrics
    sensitivity <- ifelse(TP + FN == 0, NA, TP / (TP + FN))
    specificity <- ifelse(TN + FP == 0, NA, TN / (TN + FP))
    precision   <- ifelse(TP + FP == 0, NA, TP / (TP + FP))
    f1          <- ifelse(is.na(precision) | is.na(sensitivity) | (precision + sensitivity) == 0,
                          NA, 2 * precision * sensitivity / (precision + sensitivity))
    denom_mcc   <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc         <- ifelse(denom_mcc == 0, NA, (TP * TN - FP * FN) / denom_mcc)
    jaccard     <- ifelse(TP + FP + FN == 0, NA, TP / (TP + FP + FN))
    
    # return all metrics 
    list(
      sensitivity = sensitivity, 
      specificity = specificity, 
      f1 = f1, 
      mcc = mcc, 
      jaccard = jaccard,
      n_true_edges = n_true_edges,
      n_estimated_edges = n_estimated_edges,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN
    )
  }
)



library(dplyr)
library(tidyr)
library(DT)

results_sims <- read.table("sims_batch_29.txt", header = TRUE)

# 1) Define the metrics/edge columns we care about
metric_cols <- intersect(c("sensitivity", "specificity", "f1", "mcc", "jaccard",
                           "n_true_edges", "n_estimated_edges"), names(results_sims))

# 2) Pivot longer for plotting / summarisation
results_sims_long <- results_sims %>%
  pivot_longer(all_of(metric_cols), names_to = "metric", values_to = "value")


# 3) Summarise metrics
results_sims_summmary <- results_sims_long %>%
  group_by(n_people, n_nodes, missing, estimator, metric) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")