# ------------------------------------------------------------------------------
#
# all_variable_selection.R
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


# Define all_variable_selection output folder ----------------------------------
print("Creating output/all_variable_selection output folder")

all_variable_selection_dir <- "output/all_variable_selection/"
fs::dir_create(here::here(all_variable_selection_dir))


# Load imputation object -------------------------------------------------------
print("Load imputation object")

imp <- readRDS(
  paste0("output/apply_within_MI/apply_within_MI_imp_datasets_", name, ".rds")
)


# Applying lasso, lasso_X and lasso_union models to the MI datasets
print("Applying lasso, lasso_X and lasso_union models to the MI datasets")

all_var_names <- c(
  "cov_bin_covid",
  "cov_num_age", "cov_cat_sex", "cov_num_bmi", "cov_cat_ethnicity", "cov_cat_imd", "cov_cat_smoking",
  "cov_bin_carehome", "cov_bin_hcworker", "cov_bin_dementia", "cov_bin_liver_disease", "cov_bin_ckd",
  "cov_bin_cancer", "cov_bin_hypertension", "cov_bin_diabetes", "cov_bin_obesity", "cov_bin_copd",
  "cov_bin_depression", "cov_bin_stroke_all", "cov_bin_other_ae", "cov_bin_vte", "cov_bin_hf",
  "cov_bin_angina", "cov_bin_lipidmed", "cov_bin_antiplatelet", "cov_bin_anticoagulant", "cov_bin_cocp",
  "cov_bin_hrt", "strat_cat_region"
)

lasso_var_selection_results       <- data.frame(matrix(
  nrow = get_number_of_imputed_datasets(),
  ncol = length(all_var_names)
))
rownames(lasso_var_selection_results) <- c(1:get_number_of_imputed_datasets())
colnames(lasso_var_selection_results) <- all_var_names

lasso_X_var_selection_results     <- data.frame(matrix(
  nrow = get_number_of_imputed_datasets(),
  ncol = length(all_var_names)
))
rownames(lasso_X_var_selection_results) <- c(1:get_number_of_imputed_datasets())
colnames(lasso_X_var_selection_results) <- all_var_names

lasso_union_var_selection_results <- data.frame(matrix(
  nrow = get_number_of_imputed_datasets(),
  ncol = length(all_var_names)
))
rownames(lasso_union_var_selection_results) <- c(1:get_number_of_imputed_datasets())
colnames(lasso_union_var_selection_results) <- all_var_names

for (i in c(1:get_number_of_imputed_datasets())) {
  imp_data_i <- complete(imp, action = i)

  # lasso

  lasso_cox_conf_matrix_i <- (imp_data_i %>% select(c(
    cov_bin_covid,
    cov_num_age, cov_cat_sex, cov_num_bmi, cov_cat_ethnicity, cov_cat_imd, cov_cat_smoking,
    cov_bin_carehome, cov_bin_hcworker, cov_bin_dementia, cov_bin_liver_disease, cov_bin_ckd,
    cov_bin_cancer, cov_bin_hypertension, cov_bin_diabetes, cov_bin_obesity, cov_bin_copd,
    cov_bin_depression, cov_bin_stroke_all, cov_bin_other_ae, cov_bin_vte, cov_bin_hf,
    cov_bin_angina, cov_bin_lipidmed, cov_bin_antiplatelet, cov_bin_anticoagulant, cov_bin_cocp,
    cov_bin_hrt, strat_cat_region
  )))

  lasso_cox_conf_matrix_preserving_factors_i <- model.matrix(
    as.formula(" ~ ."), # formula meaning take all terms
    data = lasso_cox_conf_matrix_i
  )
  lasso_cox_outcome_survival_i <- Surv(
    time  = as.numeric(imp_data_i$outcome_cox_dates),
    event = imp_data_i$cens_status,
    type  = "right"
  )

  cv_lasso_cox_model_i <- cv.glmnet(
    x       = lasso_cox_conf_matrix_preserving_factors_i,
    y       = lasso_cox_outcome_survival_i,
    nlambda = 200,   # length of lambda sequence
    nfolds  = get_number_of_imputed_datasets(),
    family  = "cox",      # cox regression
    alpha   = 1           # LASSO penalty
  )

  # tune regularisation parameter lambda to minimise cross-validated error (cvm)
  lambda_i <- cv_lasso_cox_model_i$lambda.min

  lasso_cox_coefs_i           <- coef(cv_lasso_cox_model_i, s = lambda_i)
  lasso_cox_coefs_i           <- as.data.frame(as.matrix(lasso_cox_coefs_i))
  colnames(lasso_cox_coefs_i) <- c("coefficient")

  non_zero_coefs_i <- lasso_cox_coefs_i %>% dplyr::filter(coefficient != 0.0)
  non_zero_vars_i  <- rownames(non_zero_coefs_i)
  lasso_vars_selected_i  <- convert_terms_to_vars(non_zero_vars_i)

  # always include exposure
  if (!("cov_bin_covid" %in% lasso_vars_selected_i)) {
    lasso_vars_selected_i <- c(lasso_vars_selected_i, "cov_bin_covid")
  }


  # lasso_X

  lasso_X_conf_matrix_i <- (imp_data_i %>% select(c(
    cov_num_age, cov_cat_sex, cov_num_bmi, cov_cat_ethnicity, cov_cat_imd, cov_cat_smoking,
    cov_bin_carehome, cov_bin_hcworker, cov_bin_dementia, cov_bin_liver_disease, cov_bin_ckd,
    cov_bin_cancer, cov_bin_hypertension, cov_bin_diabetes, cov_bin_obesity, cov_bin_copd,
    cov_bin_depression, cov_bin_stroke_all, cov_bin_other_ae, cov_bin_vte, cov_bin_hf,
    cov_bin_angina, cov_bin_lipidmed, cov_bin_antiplatelet, cov_bin_anticoagulant, cov_bin_cocp,
    cov_bin_hrt, strat_cat_region
  )))

  lasso_X_conf_matrix_preserving_factors_i <- model.matrix(
    ~ ., # formula meaning take all terms
    data = lasso_X_conf_matrix_i
  )

  lasso_X_exposure_matrix_i <- (imp_data_i %>% select(c(
    cov_bin_covid,
  )))

  lasso_X_exposure_matrix_preserving_factors_i <- model.matrix(
    ~ .,
    data = lasso_X_exposure_matrix_i
  )

  # Fitting the lasso_X logistic model ------------------------------------------
  message("Fitting the lasso_X logistic model")

  cv_lasso_X_logistic_model_i <- cv.glmnet(
    x       = lasso_X_conf_matrix_preserving_factors_i,
    y       = lasso_X_exposure_matrix_preserving_factors_i,
    nlambda = 200,        # length of lambda sequence
    nfolds  = get_number_of_imputed_datasets(),
    family  = "binomial", # logistic regression
    alpha   = 1
  )

  # tune regularisation parameter lambda to minimise cross-validated error (cvm)
  lambda_i         <- cv_lasso_X_logistic_model_i$lambda.min

  lasso_X_logistic_coefs_i           <- coef(cv_lasso_X_logistic_model_i, s = lambda_i)
  lasso_X_logistic_coefs_i           <- as.data.frame(as.matrix(lasso_X_logistic_coefs_i))
  colnames(lasso_X_logistic_coefs_i) <- c("coefficient")

  non_zero_coefs_i         <- lasso_X_logistic_coefs_i %>% dplyr::filter(coefficient != 0.0)
  non_zero_vars_i          <- rownames(non_zero_coefs_i)
  lasso_X_vars_selected_i  <- convert_terms_to_vars(non_zero_vars_i)

  # always include exposure
  if (!("cov_bin_covid" %in% lasso_X_vars_selected_i)) {
    lasso_X_vars_selected_i <- c(lasso_X_vars_selected_i, "cov_bin_covid")
  }

  # lasso_union

  lasso_union_vars_selected_i <- unique(c(lasso_vars_selected_i, lasso_X_vars_selected_i))

  # record results

  lasso_var_selection_results[i, ] <- convert_vars_to_binary_vector(
    vars   = lasso_vars_selected_i,
    labels = all_var_names
  )

  lasso_X_var_selection_results[i, ] <- convert_vars_to_binary_vector(
    vars   = lasso_X_vars_selected_i,
    labels = all_var_names
  )

  lasso_union_var_selection_results[i, ] <- convert_vars_to_binary_vector(
    vars   = lasso_union_vars_selected_i,
    labels = all_var_names
  )
}


# Aggregate variable selection results -----------------------------------------
print("Aggregate variable selection results")

lasso_mean_var_selection_results       <- colMeans(lasso_var_selection_results)
lasso_X_mean_var_selection_results     <- colMeans(lasso_X_var_selection_results)
lasso_union_mean_var_selection_results <- colMeans(lasso_union_var_selection_results)

lasso_aggregate_var_selection_results       <- names(lasso_mean_var_selection_results[lasso_mean_var_selection_results >= 0.5])
lasso_X_aggregate_var_selection_results     <- names(lasso_X_mean_var_selection_results[lasso_X_mean_var_selection_results >= 0.5])
lasso_union_aggregate_var_selection_results <- names(lasso_union_mean_var_selection_results[lasso_union_mean_var_selection_results >= 0.5])


# Save results -----------------------------------------------------------------
print("Save results")

write.csv(
  lasso_mean_var_selection_results,
  paste0(all_variable_selection_dir, "lasso_mean_var_selection_results-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  lasso_X_mean_var_selection_results,
  paste0(all_variable_selection_dir, "lasso_X_mean_var_selection_results-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  lasso_union_mean_var_selection_results,
  paste0(all_variable_selection_dir, "lasso_union_mean_var_selection_results-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  lasso_aggregate_var_selection_results,
  paste0(all_variable_selection_dir, "lasso_aggregate_var_selection_results-", name, ".csv"),
  row.names = FALSE
)

write.csv(
  lasso_X_aggregate_var_selection_results,
  paste0(all_variable_selection_dir, "lasso_X_aggregate_var_selection_results-", name, ".csv"),
  row.names = FALSE
)

write.csv(
  lasso_union_aggregate_var_selection_results,
  paste0(all_variable_selection_dir, "lasso_union_aggregate_var_selection_results-", name, ".csv"),
  row.names = FALSE
)
