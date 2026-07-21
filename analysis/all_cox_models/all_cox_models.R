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

models <- fully_adjusted_cox_models$analyses
model <- models[[1]]
print(model)
stop("how to pool????")

# pooled_model <- pool(models) # effects parametric
# 
# pool_RR <- function(est, se, conf.level=0.95, n, k){
#   m <- length(est)
#   mean_est <- mean(est)
#   var_w <-
#     mean(se^2) # within variance
#   var_b <-
#     var(est) # between variance
#   var_T <-
#     var_w + (1 + (1/m)) * var_b # total variance
#   se_total <-
#     sqrt(var_T)
#   r <- (1 + 1 / m) * (var_b / var_w)
#   v_old <- (m - 1) * (1 + (1/r))^2
#   lambda <- (var_b + (var_b/m))/var_T
#   v_obs <- (((n-k) + 1) / ((n-k) + 3)) * (n-k) * (1-lambda)
#   v_adj <- (v_old * v_obs) / (v_old + v_obs)
#   alpha <- 1 - (1 - conf.level)/2
#   t_stats <- mean_est/se_total
#   p_val <-
#     2*pt(-abs(t_stats),df=v_adj)
#   t <- qt(alpha, v_adj)
#   ci_upper <-
#     mean_est + (t*se_total)
#   ci_lower <-
#     mean_est - (t*se_total)
#   res <-
#     round(c(mean_est, ci_lower, ci_upper, p_val), 7)
#   names(res) <- c("Estimate", "95% CI L", "95% CI U", "P-val")
#   return(res)
# }
# 
# pooled_model <- pool_RR(models)
# 
# print(pooled_model)

stop("did we succeed?")

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
