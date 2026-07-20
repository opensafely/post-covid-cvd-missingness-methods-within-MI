# ------------------------------------------------------------------------------
#
# all_cox_models.R
# 
# Arguments:
#  - name - string, defines cohort ant outcome of study
#           (cohort_prevax-main-ami)
#  - cohort - string, defines which of three opensafely cohorts to describe
#             (prevax, vax, unvax)
#
# Returns:
#  - many!
#
# Authors: Emma Tarmey
#
# ------------------------------------------------------------------------------


# Load libraries ---------------------------------------------------------------
print("Load libraries")

library(magrittr)
library(mice)
library(here)
library(dplyr)
library(glmnet)
library(fs)
library(survival)


# Source common functions ------------------------------------------------------
print("Source common functions")

source("analysis/utility.R")


# Specify arguments ------------------------------------------------------------
print("Specify arguments")

args <- commandArgs(trailingOnly = TRUE)
print(length(args))

if (length(args) == 0) {
  # default argument values
  name   <- "cohort_prevax-main-ami"
  cohort <- "prevax"

} else {
  # YAML arguments
  name   <- args[[1]]
  cohort <- args[[2]]
}


# Define all_cox_models output folder -----------------------------------------
print("Creating output/all_cox_models output folder")

all_cox_models_dir <- "output/all_cox_models/"
fs::dir_create(here::here(all_cox_models_dir))


# Load data -------------------------------------------------------------------
print("Load data")

imp <- readRDS(
  paste0("output/apply_within_MI/apply_within_MI_imp_datasets_", name, ".rds")
)

fully_adjusted_var_selection_results <- c(
  "cov_bin_ami",
  "cov_bin_sahhs",
  # "cov_bin_covid", # avoid double counting

  # "cov_num_age", # avoid double counting
  # "cov_cat_sex", # avoid double counting
  "cov_num_bmi",
  "cov_cat_ethnicity",
  "cov_cat_imd",
  "cov_cat_smoking",

  "cov_bin_carehome",
  "cov_bin_hcworker",
  "cov_bin_dementia",
  "cov_bin_liver_disease",
  "cov_bin_ckd",

  "cov_bin_cancer",
  "cov_bin_hypertension",
  "cov_bin_diabetes",
  "cov_bin_obesity",
  "cov_bin_copd",

  "cov_bin_depression",
  "cov_bin_stroke_all",
  "cov_bin_other_ae",
  "cov_bin_vte",
  "cov_bin_hf",

  "cov_bin_angina",
  "cov_bin_lipidmed",
  "cov_bin_antiplatelet",
  "cov_bin_anticoagulant",
  "cov_bin_cocp",

  "cov_bin_hrt"
  # "strat_cat_region" # avoid double counting
)

# remove outcome from variable selection
if (grepl("ami", name)) {
  fully_adjusted_var_selection_results <- fully_adjusted_var_selection_results[! fully_adjusted_var_selection_results %in% c("cov_bin_ami")]
} else{
  fully_adjusted_var_selection_results <- fully_adjusted_var_selection_results[! fully_adjusted_var_selection_results %in% c("cov_bin_sahhs")]
}

lasso_var_selection_results <- read.csv(
  paste0("output/all_variable_selection/lasso_aggregate_var_selection_results-", name, ".csv")
)[, 'x']

lasso_X_var_selection_results <- read.csv(
  paste0("output/all_variable_selection/lasso_X_aggregate_var_selection_results-", name, ".csv")
)[, 'x']

lasso_union_var_selection_results <- read.csv(
  paste0("output/all_variable_selection/lasso_union_aggregate_var_selection_results-", name, ".csv")
)[, 'x']


# do not souble count age, sex, binary exposure or region
double_counts                     <- c("cov_cat_sex", "cov_num_age", "cov_bin_covid", "strat_cat_region")
lasso_var_selection_results       <- lasso_var_selection_results[! lasso_var_selection_results %in% double_counts]
lasso_X_var_selection_results     <- lasso_X_var_selection_results[! lasso_X_var_selection_results %in% double_counts]
lasso_union_var_selection_results <- lasso_union_var_selection_results[! lasso_union_var_selection_results %in% double_counts]


# Fit cox regression models on imputed datasets -----------------------------
print("Fit cox regression models on imputed datasets")

source("analysis/cox_ipw/fn-survival_data_setup.R")
source("analysis/cox_ipw/fn-get_episode_info.R")
source("analysis/cox_ipw/fn-fit_model.R")
source("analysis/cox_ipw/fn-preprocess_and_fit_model.R")

for (i in c(1:get_number_of_imputed_datasets())) {
  model_input_df_i <- complete(imp, action = i)

  fully_adjusted_cox_model_i <- preprocess_and_fit_model(
    df              = model_input_df_i,
    covariate_other = fully_adjusted_var_selection_results
  )

  lasso_cox_model_i <- preprocess_and_fit_model(
    df              = model_input_df_i,
    covariate_other = lasso_var_selection_results
  )

  lasso_X_cox_model_i <- preprocess_and_fit_model(
    df              = model_input_df_i,
    covariate_other = lasso_X_var_selection_results
  )

  lasso_union_cox_model_i <- preprocess_and_fit_model(
    df              = model_input_df_i,
    covariate_other = lasso_union_var_selection_results
  )
}

# NB: warnings refer to constant columns, which are intended and can be safely ignored
stop("pool the above!")


# Save results -------------------------------------------------------------
print("Save results")

write.csv(
  pooled_fully_adjusted_cox_model,
  paste0(all_cox_models_dir, "pooled_fully_adjusted_cox_model-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  pooled_lasso_cox_model,
  paste0(all_cox_models_dir, "pooled_lasso_cox_model-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  pooled_lasso_X_cox_model,
  paste0(all_cox_models_dir, "pooled_lasso_X_cox_model-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  pooled_lasso_union_cox_model,
  paste0(all_cox_models_dir, "pooled_lasso_union_cox_model-", name, ".csv"),
  row.names = TRUE
)
