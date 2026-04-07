library(tidyverse)
library(dplyr)

# 1) COMBINE ALL BATCHES
all_batches <- lapply(1:100, function(i) {
  fname <- paste0("sims_batch_", i, ".txt")
  if (file.exists(fname)) {
    read.table(fname, header = TRUE)
  } else {
    message("Missing: ", fname)
    NULL
  }
}) %>% bind_rows()

# Check
all_batches %>% count(n_people, missing, estimator) %>% pull(n) %>% unique()
# Should return 200


# 2) DEFINE METRIC COLUMNS
metric_cols <- intersect(c("sensitivity", "specificity", "f1", "mcc", "jaccard",
                           "n_true_edges", "n_estimated_edges"), names(all_batches))

# 3) PIVOT LONGER
all_batches_long <- all_batches %>%
  pivot_longer(all_of(metric_cols), names_to = "metric", values_to = "value")

# 4) SUMMARISE
all_batches_summary <- all_batches_long %>%
  group_by(n_people, n_nodes, missing, estimator, metric) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    n_reps     = n(),
    .groups = "drop"
  )


view(all_batches_summary)

write.csv(all_batches_summary, "all_batches_summary.csv", row.names = FALSE)




# family comparison

all_batches_summary <- read.csv("all_batches_summary.csv")

# family groupings 
family_map <- tribble(
  ~estimator,                          ~family,
  "EBICglasso_listwise",               "Deletion-based",
  "EBICglasso_pairwise",               "Deletion-based",
  "ggmModSelect_listwise",             "Deletion-based",
  "ggmModSelect_pairwise",             "Deletion-based",
  "EBICglasso_mean",                   "Simple imputation",
  "EBICglasso_median",                 "Simple imputation",
  "EBICglasso_knn",                    "ML-based",
  "EBICglasso_single_imputation_pmm",  "Simple imputation",
  "ggmModSelect_mean",                 "Simple imputation",
  "ggmModSelect_median",               "Simple imputation",
  "ggmModSelect_knn",                  "ML-based",
  "ggmModSelect_single_imputation_pmm","Simple imputation",
  "EBICglasso_miceMI",                 "MICE",
  "ggmModSelect_miceMI",               "MICE",
  "EBICglasso_missForest",             "ML-based",
  "EBICglasso_xgboost",                "ML-based",
  "ggmModSelect_missForest",           "ML-based",
  "EBICglasso_fiml",                   "Model-based",
  "ggmModSelect_fiml",                 "Model-based",
  "ggm_prune",                         "Model-based",
  "EM_EBIC",                           "Model-based",
  "EM_CV",                             "Model-based"
)

summary_families <- all_batches_summary %>%
  left_join(family_map, by = "estimator")

# Compute family means (averaged across estimators AND sample sizes)
family_trends <- summary_families %>%
  filter(metric %in% c("sensitivity", "specificity", "mcc")) %>%
  group_by(family, missing, metric) %>%
  summarise(
    family_mean = mean(mean_value, na.rm = TRUE),
    .groups = "drop"
  )

family_trends %>%
  filter(metric == "mcc") %>%
  pivot_wider(names_from = missing, values_from = family_mean) %>%
  arrange(family) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  view()


family_trends %>%
  filter(metric == "sensitivity") %>%
  pivot_wider(names_from = missing, values_from = family_mean) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  view()

family_trends %>%
  filter(metric == "specificity") %>%
  pivot_wider(names_from = missing, values_from = family_mean) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  view()

procedure_comparison <- summary_families %>%
  filter(metric %in% c("sensitivity", "specificity", "mcc")) %>%
  mutate(procedure = case_when(
    startsWith(estimator, "EBICglasso")   ~ "EBICglasso",
    startsWith(estimator, "ggmModSelect") ~ "ggmModSelect",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(procedure)) %>%
  group_by(procedure, missing, metric) %>%
  summarise(proc_mean = mean(mean_value, na.rm = TRUE), .groups = "drop")

procedure_comparison %>%
  pivot_wider(names_from = procedure, values_from = proc_mean) %>%
  mutate(gap = round(ggmModSelect - EBICglasso, 3)) %>%
  filter(metric %in% c("specificity", "mcc")) %>%
  select(metric, missing, EBICglasso, ggmModSelect, gap) %>%
  arrange(metric, missing) %>%
  view()

summary_families %>%
  filter(metric == "sensitivity") %>%
  group_by(n_people, missing) %>%
  summarise(mean_sens = mean(mean_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = n_people, values_from = mean_sens) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  view()


composite <- all_batches_summary %>%
  filter(metric %in% c("sensitivity", "specificity", "mcc")) %>%
  group_by(estimator, metric) %>%
  summarise(overall_mean = mean(mean_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = metric, values_from = overall_mean) %>%
  mutate(
    composite = (sensitivity + specificity + mcc) / 3,
    across(where(is.numeric), ~round(., 3))
  ) %>%
  arrange(desc(composite))

view(composite)
