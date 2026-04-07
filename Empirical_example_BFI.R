library(bootnet)
library(qgraph)
library(psychonetrics)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mice)
library(EMgaussian)
library(missForest)
library(VIM)
library(mixgb)
library(psych)

# 1. LOAD DATA
data(bfi)

# 8 items: 4 Neuroticism + 4 Conscientiousness
items     <- c("N1", "N2", "N3", "N4", "C1", "C2", "C3", "C4")
bfi_sub   <- bfi[, items]
bfi_clean <- na.omit(bfi_sub)
cat("N after removing existing NAs:", nrow(bfi_clean), "\n")

# 2. TRUE NETWORK (complete data, used as ground truth)

modS_true   <- varcov(bfi_clean, estimator = "FIML") %>% runmodel()
S_true      <- getmatrix(modS_true, "sigma")
trueNetwork <- qgraph::ggmModSelect(S_true, n = nrow(bfi_clean))$graph
cat("True edges:", sum(trueNetwork[lower.tri(trueNetwork)] != 0), "\n")


# 3. IMPOSE MCAR

impose_mcar <- function(data, prop) {
  data_miss <- data
  n <- nrow(data) # number of participants e.g. 2700
  p <- ncol(data) # number of variables e.g. 8
  n_missing <- round(prop * n * p) # % of cells to delete
  idx  <- sample(n * p, n_missing) # randomy pick some cells
  rows <- ((idx - 1) %% n) + 1  # convert each flat index into a row and column position
  cols <- ((idx - 1) %/% n) + 1
  for (i in seq_along(rows)) data_miss[rows[i], cols[i]] <- NA
  return(data_miss) # Loop through each selected cell and set it to NA.
}

# 4. RECOVERY METRICS (identical to simulation code)

compute_metrics <- function(trueNetwork, estNetwork) {
  lt <- lower.tri(trueNetwork, diag = FALSE)
  tE <- trueNetwork[lt]
  eE <- estNetwork[lt]
  
  TP <- sum(tE != 0 & eE != 0)
  FP <- sum(tE == 0 & eE != 0)
  FN <- sum(tE != 0 & eE == 0)
  TN <- sum(tE == 0 & eE == 0)
  
  sensitivity <- ifelse(TP + FN == 0, NA, TP / (TP + FN))
  specificity <- ifelse(TN + FP == 0, NA, TN / (TN + FP))
  precision   <- ifelse(TP + FP == 0, NA, TP / (TP + FP))
  f1          <- ifelse(is.na(precision) | is.na(sensitivity) | (precision + sensitivity) == 0,
                        NA, 2 * precision * sensitivity / (precision + sensitivity))
  denom_mcc   <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc         <- ifelse(denom_mcc == 0, NA, (TP * TN - FP * FN) / denom_mcc)
  jaccard     <- ifelse(TP + FP + FN == 0, NA, TP / (TP + FP + FN))
  
  list(sensitivity = sensitivity, specificity = specificity,
       f1 = f1, mcc = mcc, jaccard = jaccard,
       n_true_edges      = sum(tE != 0),
       n_estimated_edges = sum(eE != 0))
}

# 5. RUN ALL 22 ESTIMATORS

miss_levels <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
n_reps      <- 5 # just 5 reps
threshold   <- 3  # same as simulation

all_estimators <- c(
  "ggm_prune",
  "EBICglasso_listwise",
  "EBICglasso_pairwise",
  "EBICglasso_fiml",
  "EBICglasso_single_imputation_pmm",
  "EBICglasso_miceMI",
  "ggmModSelect_listwise",
  "ggmModSelect_pairwise",
  "ggmModSelect_fiml",
  "ggmModSelect_single_imputation_pmm",
  "ggmModSelect_miceMI",
  "EM_EBIC",
  "EM_CV",
  "EBICglasso_mean",
  "EBICglasso_median",
  "EBICglasso_knn",
  "EBICglasso_missForest",
  "ggmModSelect_mean",
  "ggmModSelect_median",
  "ggmModSelect_knn",
  "ggmModSelect_missForest",
  "EBICglasso_xgboost"
)

results_list <- list()
set.seed(67)

for (missing in miss_levels) {
  cat("\n── Missingness:", missing, "──\n")
  
  for (rep in 1:n_reps) {
    
    simData <- impose_mcar(bfi_clean, missing)
    
    for (estimator in all_estimators) {
      
      estNetwork <- tryCatch({
        
        if (estimator == "ggm_prune") {
          mod <- ggm(simData, estimator = "FIML") %>% runmodel()
          SEs <- !any(is.na(mod@parameters$se[mod@parameters$matrix == "omega"]))
          if (SEs) {
            mod <- mod %>% prune(alpha = 0.05)
            getmatrix(mod, "omega")
          } else stop("SEs could not be computed")
          
        } else if (estimator == "EBICglasso_listwise") {
          estimateNetwork(simData, default = "EBICglasso", missing = "listwise")$graph
          
        } else if (estimator == "EBICglasso_pairwise") {
          estimateNetwork(simData, default = "EBICglasso", missing = "pairwise")$graph
          
        } else if (estimator == "EBICglasso_fiml") {
          modS <- varcov(simData, estimator = "FIML") %>% runmodel()
          S_fiml <- getmatrix(modS, "sigma")
          qgraph::EBICglasso(S_fiml, n = nrow(simData))
          
        } else if (estimator == "EBICglasso_single_imputation_pmm") {
          imp <- mice(simData, m = 1, method = "pmm", maxit = 5, printFlag = FALSE)
          dat_mi <- complete(imp, 1)
          estimateNetwork(dat_mi, default = "EBICglasso")$graph
          
        } else if (estimator == "EBICglasso_miceMI") {
          imp <- mice(simData, m = 5, method = "pmm", maxit = 5, printFlag = FALSE)
          graphs <- lapply(1:5, function(i) {
            dat_i <- complete(imp, i)
            estimateNetwork(dat_i, default = "EBICglasso")$graph
          })
          arr <- simplify2array(graphs)
          apply(arr, c(1, 2), function(x) {
            nonzero <- x[x != 0]
            if (length(nonzero) >= threshold) median(nonzero) else 0
          })
          
        } else if (estimator == "ggmModSelect_listwise") {
          estimateNetwork(simData, default = "ggmModSelect", missing = "listwise")$graph
          
        } else if (estimator == "ggmModSelect_pairwise") {
          estimateNetwork(simData, default = "ggmModSelect", missing = "pairwise")$graph
          
        } else if (estimator == "ggmModSelect_fiml") {
          modS <- varcov(simData, estimator = "FIML") %>% runmodel()
          S_fiml <- getmatrix(modS, "sigma")
          qgraph::ggmModSelect(S_fiml, n = nrow(simData))$graph
          
        } else if (estimator == "ggmModSelect_single_imputation_pmm") {
          imp <- mice(simData, m = 1, method = "pmm", maxit = 5, printFlag = FALSE)
          dat_mi <- complete(imp, 1)
          estimateNetwork(dat_mi, default = "ggmModSelect")$graph
          
        } else if (estimator == "ggmModSelect_miceMI") {
          imp <- mice(simData, m = 5, method = "pmm", maxit = 5, printFlag = FALSE)
          graphs <- lapply(1:5, function(i) {
            dat_i <- complete(imp, i)
            estimateNetwork(dat_i, default = "ggmModSelect")$graph
          })
          arr <- simplify2array(graphs)
          apply(arr, c(1, 2), function(x) {
            nonzero <- x[x != 0]
            if (length(nonzero) >= threshold) median(nonzero) else 0
          })
          
        } else if (estimator == "EM_EBIC") {
          X   <- scale(as.matrix(simData))
          rho <- EMgaussian::rhogrid(100, method = "qgraph", dat = simData)
          estimateNetwork(X, fun = EMgaussian::EMggm, rho = rho,
                          glassoversion = "glasso", rhoselect = "ebic")$graph
          
        } else if (estimator == "EM_CV") {
          X   <- scale(as.matrix(simData))
          rho <- EMgaussian::rhogrid(100, method = "qgraph", dat = X)
          estimateNetwork(X, fun = EMgaussian::EMggm, rho = rho,
                          glassoversion = "glasso", rhoselect = "kfold", k = 5)$graph
          
        } else if (estimator == "EBICglasso_mean") {
          dat_imp <- as.data.frame(simData)
          for (j in seq_along(dat_imp)) {
            if (is.numeric(dat_imp[[j]])) {
              na_idx <- is.na(dat_imp[[j]])
              dat_imp[[j]][na_idx] <- mean(dat_imp[[j]], na.rm = TRUE)
            }
          }
          estimateNetwork(dat_imp, default = "EBICglasso")$graph
          
        } else if (estimator == "EBICglasso_median") {
          dat_imp <- as.data.frame(simData)
          for (j in seq_along(dat_imp)) {
            if (is.numeric(dat_imp[[j]])) {
              na_idx <- is.na(dat_imp[[j]])
              dat_imp[[j]][na_idx] <- stats::median(dat_imp[[j]], na.rm = TRUE)
            }
          }
          estimateNetwork(dat_imp, default = "EBICglasso")$graph
          
        } else if (estimator == "EBICglasso_knn") {
          dat_imp <- VIM::kNN(simData, k = 5, imp_var = FALSE)
          estimateNetwork(dat_imp, default = "EBICglasso")$graph
          
        } else if (estimator == "EBICglasso_missForest") {
          dat_imp <- missForest(as.data.frame(simData), maxiter = 5, ntree = 100, verbose = FALSE)$ximp
          estimateNetwork(dat_imp, default = "EBICglasso")$graph
          
        } else if (estimator == "ggmModSelect_mean") {
          dat_imp <- as.data.frame(simData)
          for (j in seq_along(dat_imp)) {
            if (is.numeric(dat_imp[[j]])) {
              na_idx <- is.na(dat_imp[[j]])
              dat_imp[[j]][na_idx] <- mean(dat_imp[[j]], na.rm = TRUE)
            }
          }
          estimateNetwork(dat_imp, default = "ggmModSelect")$graph
          
        } else if (estimator == "ggmModSelect_median") {
          dat_imp <- as.data.frame(simData)
          for (j in seq_along(dat_imp)) {
            if (is.numeric(dat_imp[[j]])) {
              na_idx <- is.na(dat_imp[[j]])
              dat_imp[[j]][na_idx] <- stats::median(dat_imp[[j]], na.rm = TRUE)
            }
          }
          estimateNetwork(dat_imp, default = "ggmModSelect")$graph
          
        } else if (estimator == "ggmModSelect_knn") {
          dat_imp <- VIM::kNN(simData, k = 5, imp_var = FALSE)
          estimateNetwork(dat_imp, default = "ggmModSelect")$graph
          
        } else if (estimator == "ggmModSelect_missForest") {
          dat_imp <- missForest(as.data.frame(simData), maxiter = 5, ntree = 100, verbose = FALSE)$ximp
          estimateNetwork(dat_imp, default = "ggmModSelect")$graph
          
        } else if (estimator == "EBICglasso_xgboost") {
          dat_xgb <- as.data.frame(simData)
          dat_xgb[] <- lapply(dat_xgb, as.numeric)
          dat_imp <- mixgb(data = dat_xgb, m = 1, nrounds = 50, maxit = 1, verbose = FALSE,
                           xgb.params = list(max_depth = 4, eta = 0.3, nthread = 1))
          dat_imp <- as.data.frame(dat_imp)
          estimateNetwork(dat_imp, default = "EBICglasso")$graph
          
        } else {
          stop("Invalid estimator")
        }
        
      }, error = function(e) {
        message("  Error [", estimator, " | miss=", missing, " | rep=", rep, "]: ", e$message)
        NULL
      })
      
      if (!is.null(estNetwork)) {
        m <- compute_metrics(trueNetwork, estNetwork)
        results_list[[length(results_list) + 1]] <- data.frame(
          missing           = missing,
          rep               = rep,
          estimator         = estimator,
          sensitivity       = m$sensitivity,
          specificity       = m$specificity,
          f1                = m$f1,
          mcc               = m$mcc,
          jaccard           = m$jaccard,
          n_true_edges      = m$n_true_edges,
          n_estimated_edges = m$n_estimated_edges
        )
      }
    }
  }
  cat("  Missingness", missing, "complete\n")
}

# 6. SUMMARISE (same structure as all_batches_summary.csv)

results <- bind_rows(results_list)

metric_cols <- c("sensitivity", "specificity", "f1", "mcc", "jaccard",
                 "n_true_edges", "n_estimated_edges")

empirical_summary <- results %>%
  pivot_longer(all_of(metric_cols), names_to = "metric", values_to = "value") %>%
  group_by(estimator, missing, metric) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value,  na.rm = TRUE),
    n_reps     = n(),
    .groups    = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

print(empirical_summary)
View(empirical_summary)
write.csv(empirical_summary, "empirical_summary.csv", row.names = FALSE)