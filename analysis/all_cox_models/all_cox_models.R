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
  "cov_bin_covid",

  "cov_num_age",
  "cov_cat_sex",
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

  "cov_bin_hrt",
  "strat_cat_region"
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


# always include age, sex
lasso_var_selection_results       <- union(lasso_var_selection_results,       c("cov_cat_sex", "cov_num_age"))
lasso_X_var_selection_results     <- union(lasso_X_var_selection_results,     c("cov_cat_sex", "cov_num_age"))
lasso_union_var_selection_results <- union(lasso_union_var_selection_results, c("cov_cat_sex", "cov_num_age"))

fully_adjusted_model_formula <- make_outcome_formula(vars_selected = fully_adjusted_var_selection_results)
lasso_model_formula          <- make_outcome_formula(vars_selected = lasso_var_selection_results)
lasso_X_model_formula        <- make_outcome_formula(vars_selected = lasso_X_var_selection_results)
lasso_union_model_formula    <- make_outcome_formula(vars_selected = lasso_union_var_selection_results)


# Fit cox regression models on imputed datasets -----------------------------
print("Fit cox regression models on imputed datasets")

# TODO: ALL PREPROCESSING FROM COX-IPW.R HERE
# ABANDON WITH() PARADIGM UNLESS VECTORISES LATER

source("analysis/cox_ipw/fn-survival_data_setup.R")
source("analysis/cox_ipw/fn-get_episode_info.R")
source("analysis/cox_ipw/fn-fit_model.R")

model_input_df             <- complete(imp, action = 1)
model_input_df$study_start <- NULL
model_input_df$study_stop  <- NULL

cut_points         <- c(1, 28, 196, 364, 714, 1582)
time_period_labels <- c(
  "days0_1", "days1_28", "days28_196", "days196_364", "days364_714",
  "days714_1582"
)

episode_labels <- data.frame(
  episode          = 0:length(cut_points),
  time_period      = c("days_pre", time_period_labels),
  stringsAsFactors = FALSE
)

data_surv <- survival_data_setup(
  df             = model_input_df,
  cut_points     = cut_points,
  episode_labels = episode_labels
)

stop("got here")

episode_info <- get_episode_info(
  df             = data_surv,
  cut_points     = cut_points,
  episode_labels = episode_labels,
  ipw            = FALSE
)

example_model <- fit_model(
  df                  = data_surv,
  time_periods        = episode_info[
    episode_info$time_period != "days_pre",
  ]$time_period,
  covariates          = covariate_other,
  strata              = "strat_cat_region",
  age_spline          = TRUE,
  covariate_removed   = NULL,
  covariate_collapsed = NULL,
  ipw                 = FALSE
)

stop("above????")

# fully_adjusted
fully_adjusted_cox_models <- with(
  data = imp,
  exp  = coxph(formula = as.formula(fully_adjusted_model_formula))
)
pooled_fully_adjusted_cox_model <- summary(pool(fully_adjusted_cox_models))

# lasso
lasso_cox_models <- with(
  data = imp,
  exp  = coxph(formula = as.formula(lasso_model_formula))
)
pooled_lasso_cox_model <- summary(pool(lasso_cox_models))

# lasso_X
lasso_X_cox_models <- with(
  data = imp,
  exp  = coxph(formula = as.formula(lasso_X_model_formula))
)
pooled_lasso_X_cox_model <- summary(pool(lasso_X_cox_models))

# lasso_union
lasso_union_cox_models <- with(
  data = imp,
  exp  = coxph(formula = as.formula(lasso_union_model_formula))
)
pooled_lasso_union_cox_model <- summary(pool(lasso_union_cox_models))


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
