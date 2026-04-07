# PL4401-Thesis-Xavier
PL4401 Honours Thesis R Scripts
**Author:** Tan Xuan Rong Xavier | **Supervisor:** Dr Sacha Epskamp | NUS Psychology AY2025/2026

# 1. Simulation-HT

Monte Carlo simulation study for PL4401 Honours Thesis:
*Evaluating Imputation Strategies for Missing Data in Network Analysis: A Simulation-Based Approach*

## What this script does
Runs a full-factorial simulation comparing 22 combinations of missing-data 
handling strategies and GGM estimation procedures across:
- 3 sample sizes: N = 250, 500, 750
- 10 missingness levels: 0% to 90% (MCAR)
- 200 replications per condition (660 unique conditions total)

For each replication, a true 8-node GGM is generated (edge density = 0.2, 
80% positive edges), data are drawn from it, missingness is imposed, and 
each estimator attempts to recover the true network.

## Output
Writes batch result files (`sims_batch_X.txt`) with one row per replication, 
containing sensitivity, specificity, F1, MCC, Jaccard, edge counts, TP, FP, FN, TN.

## Requirements
```r
install.packages(c("bootnet", "qgraph", "psychonetrics", "parSim",
                   "mice", "EMgaussian", "missForest", "VIM", "mixgb",
                   "dplyr", "tidyr", "ggplot2", "tidyverse", "Matrix"))
```

## How to run
Run the script directly in R. Adjust `nCores` to match your machine.
Results are written to your working directory as `sims_batch_1.txt` through 
`sims_batch_N.txt` depending on how you batch the runs.

------------------------------------------------------------------------------------------------------------

# 2. Empirical-Example-BFI-HT

Empirical illustration for PL4401 Honours Thesis applying all 22 estimators 
to real Big Five Inventory (BFI) data to validate simulation findings.

## What this script does
Applies all 22 network estimation procedures to 8 BFI items (N1-N4, C1-C4) 
from the psych R package (N = 2,649 after listwise deletion of existing NAs).

## Output
- `empirical_summary.csv` — mean sensitivity, specificity, F1, MCC, Jaccard 
  per estimator and missingness level (averaged across 5 replications)

## Requirements
```r
install.packages(c("bootnet", "qgraph", "psychonetrics", "mice",
                   "EMgaussian", "missForest", "VIM", "mixgb", "psych",
                   "dplyr", "tidyr", "ggplot2"))
```

## How to run
Run the full script in R. Output CSV is saved to your working directory.
Runtime is approximately 30-60 minutes depending on your machine.

------------------------------------------------------------------------------------------------------------

# 3. Results-HT

Data processing and results summarisation script 

## What this script does
1. Combines all simulation batch files (`sims_batch_1.txt` to `sims_batch_100.txt`) 
   into a single data frame
2. Pivots and summarises results (mean, SD) per estimator, sample size, 
   missingness level, and metric
3. Computes family-level comparisons (deletion, simple imputation, MICE, 
   ML-based, model-based)
4. Computes EBICglasso vs ggmModSelect procedure-level comparisons
5. Generates composite performance rankings across estimators

## Output
- `all_batches_summary.csv` — primary summary file used by the Shiny dashboard
- Various summary tables printed to console and opened in the RStudio viewer

## Requirements
```r
install.packages(c("dplyr", "tidyr", "tidyverse", "DT"))
```
## How to run
Ensure all `sims_batch_X.txt` files are in your working directory, 
then run the script. The final summary CSV is saved to your working directory 
and can be loaded directly into the Shiny dashboard.

------------------------------------------------------------------------------------------------------------

# 4. Thesis-Dashboard-HT

Interactive R Shiny dashboard for PL4401 Honours Thesis simulation results.

**Live version:** https://xaviertanxr.shinyapps.io/PL4401-Thesis-Dashboard-Xavier/

## What this dashboard does
Visualises all simulation results across 22 estimators, 3 sample sizes, 
and 10 missingness levels (0-90%). 

## Requirements
```r
install.packages(c("shiny", "ggplot2", "plotly", "dplyr", "tidyr", "DT"))
```
------------------------------------------------------------------------------------------------------------
