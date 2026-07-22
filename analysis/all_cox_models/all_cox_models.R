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
source("analysis/all_cox_models/fn-get_cox_results.R")
source("analysis/all_cox_models/fn-manual_pool_RR.R")


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

imp_surv <- readRDS(
  paste0("output/apply_within_MI/apply_within_MI_imp_surv_datasets_", name, ".rds")
)

episode_info <- readRDS(
  paste0("output/apply_within_MI/apply_within_MI_episode_info_", name, ".rds")
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


fully_adjusted_formula <- make_rms_cph_formula(
  vars_selected = fully_adjusted_var_selection_results,
  time_periods  = episode_info[episode_info$time_period != "days_pre",]$time_period,
  strata        = NULL
)

lasso_formula <- make_rms_cph_formula(
  vars_selected = lasso_var_selection_results,
  time_periods  = episode_info[episode_info$time_period != "days_pre",]$time_period,
  strata        = NULL
)

lasso_X_formula <- make_rms_cph_formula(
  vars_selected = lasso_X_var_selection_results,
  time_periods  = episode_info[episode_info$time_period != "days_pre",]$time_period,
  strata        = NULL
)

lasso_union_formula <- make_rms_cph_formula(
  vars_selected = lasso_union_var_selection_results,
  time_periods  = episode_info[episode_info$time_period != "days_pre",]$time_period,
  strata        = NULL
)


# Fit cox regression models on imputed datasets -----------------------------
print("Fit cox regression models on imputed datasets")

model_input_surv_data_1 <- complete(imp_surv, action = 1)
dd <<- rms::datadist(model_input_surv_data_1)

withr::local_options(list(
  datadist = "dd",
  contrasts = c("contr.treatment", "contr.treatment")
))

fully_adjusted_cox_models <- with(
  data = imp_surv,
  expr = rms::cph(
    formula = as.formula(fully_adjusted_formula),
    # data    = environment(data),
    method  = "breslow",
    surv    = TRUE,
    x       = FALSE,
    y       = FALSE
  )
)

lasso_cox_models <- with(
  data = imp_surv,
  expr = rms::cph(
    formula = as.formula(lasso_formula),
    # data    = environment(data),
    method  = "breslow",
    surv    = TRUE,
    x       = FALSE,
    y       = FALSE
  )
)

lasso_X_cox_models <- with(
  data = imp_surv,
  expr = rms::cph(
    formula = as.formula(lasso_X_formula),
    # data    = environment(data),
    method  = "breslow",
    surv    = TRUE,
    x       = FALSE,
    y       = FALSE
  )
)

lasso_union_cox_models <- with(
  data = imp_surv,
  expr = rms::cph(
    formula = as.formula(lasso_union_formula),
    # data    = environment(data),
    method  = "breslow",
    surv    = TRUE,
    x       = FALSE,
    y       = FALSE
  )
)

fully_adjusted_cox_results <- get_cox_results(list_of_cox_models = fully_adjusted_cox_models$analyses, surv_formula = fully_adjusted_formula)
lasso_cox_results          <- get_cox_results(list_of_cox_models = lasso_cox_models$analyses,          surv_formula = lasso_formula)
lasso_X_cox_results        <- get_cox_results(list_of_cox_models = lasso_X_cox_models$analyses,        surv_formula = lasso_X_formula)
lasso_union_cox_results    <- get_cox_results(list_of_cox_models = lasso_union_cox_models$analyses,    surv_formula = lasso_union_formula)

fully_adjusted_pooled_cox_results <- manual_pool_RR(list_of_cox_results = fully_adjusted_cox_results)
lasso_pooled_cox_results          <- manual_pool_RR(list_of_cox_results = lasso_cox_results)
lasso_X_pooled_cox_results        <- manual_pool_RR(list_of_cox_results = lasso_X_cox_results)
lasso_union_pooled_cox_results    <- manual_pool_RR(list_of_cox_results = lasso_union_cox_results)


# Save results -------------------------------------------------------------
print("Save results")

write.csv(
  fully_adjusted_pooled_cox_results,
  paste0(all_cox_models_dir, "fully_adjusted_pooled_cox_results-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  lasso_pooled_cox_results,
  paste0(all_cox_models_dir, "lasso_pooled_cox_results-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  lasso_X_pooled_cox_results,
  paste0(all_cox_models_dir, "lasso_X_pooled_cox_results-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  lasso_union_pooled_cox_results,
  paste0(all_cox_models_dir, "lasso_union_pooled_cox_results-", name, ".csv"),
  row.names = TRUE
)
