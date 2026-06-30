# ------------------------------------------------------------------------------
#
# apply_within_MI.R
#
# This file applies multiple imputation to the BMI and Smoking covariates
# MI is conducted in "within" fashion, meaning that the subsequent analysis
# models are fit in parallel across all imputed datasets, results then
# being aggregated across these parallel models
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


# Define apply_within_MI output folder -------------------------------------
print("Creating output/apply_within_MI output folder")

apply_within_MI_dir <- "output/apply_within_MI/"
fs::dir_create(here::here(apply_within_MI_dir))


# Load data ----------------------------------------------------------
print("Load data")

df <- readr::read_rds(paste0(
  "output/dataset_clean/input_",
  cohort,
  "_clean.rds"
))

model_input_df <- readr::read_rds(paste0(
  "output/model/model_input-",
  name,
  ".rds"
))


# Applying multiple imputation to BMI and smoking covariates -------------------
print("Applying multiple imputation to BMI and smoking covariates")

# set random seed
set.seed(2026)

# re-cast missing smoking to NA
smoking_missing <- model_input_df$cov_cat_smoking == "Missing"
model_input_df$cov_cat_smoking[smoking_missing] <- NA

# check missingness of smoking and bmi variables
percent_smoking_missing <- signif(100 * (sum(is.na(model_input_df$cov_cat_smoking)) / length(model_input_df$cov_cat_smoking)), digits = 4)
percent_bmi_missing     <- signif(100 * (sum(is.na(model_input_df$cov_num_bmi))     / length(model_input_df$cov_num_bmi)),     digits = 4)

print(paste0("The variable smoking is ", percent_smoking_missing, "% missing"))
print(paste0("The variable bmi is ",     percent_bmi_missing,     "% missing"))


# Specify imputation methods for each outcome (ami and sahhs) ------------------
print("Specify imputation methods for each outcome (ami and sahhs)")

# The below ensures that bmi and smoking are handled with specific
# imputation methods and that all other covariates are left alone
# See: https://www.rdocumentation.org/packages/mice/versions/3.17.0/topics/mice
imp_method                    <- rep("", length.out = length(colnames(model_input_df)))
names(imp_method)             <- colnames(model_input_df)
imp_method["cov_cat_smoking"] <- "polyreg" # smoking is categorical with 3 levels, Polytomous logistic regression
imp_method["cov_num_bmi"]     <- "norm"    # bmi is numerical, Bayesian linear regression


# Specify imputation formulas for each outcome (ami and sahhs) -----------------
print("Specify imputation formulas for outcome")

# Specify imputation formulas, exclude variable such as index date
all_var_names <- c(
  "cov_bin_ami",
  "cov_bin_sahhs",
  "cov_bin_covid",

  "cov_num_age",
  "cov_cat_sex",
  # "cov_num_bmi", # excluded
  "cov_cat_ethnicity",
  "cov_cat_imd",
  # "cov_cat_smoking", # excluded

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
  # "vax_cat_jcvi_group" # excluded
  # "cens_status" # excluded
)

my_formulas <- list(
  cov_cat_smoking = as.formula(paste0("cov_cat_smoking ~ ", paste(all_var_names, collapse = " + "), " + H0")),
  cov_num_bmi     = as.formula(paste0("cov_num_bmi ~ ",     paste(all_var_names, collapse = " + "), " + H0"))
)


# Calculate Nelson-Aalen Estimator for outcome -----------
print("Calculate Nelson-Aalen Estimator for outcome")

if (grepl("ami", name)) {
  outcome <- "ami"
} else {
  outcome <- "sahhs"
}

outcome_cox_dates <- rep(as.Date(NA), times = nrow(model_input_df))
cens_status       <- rep(NA, times = nrow(model_input_df))

# 0 = censoring time = date of end of study
# 1 = failure time = time of outcome event
# See: https://www.rdocumentation.org/packages/survival/versions/3.8-3/topics/Surv
# and: https://glmnet.stanford.edu/articles/Coxnet.html#basic-usage-for-right-censored-data
for (i in c(1:nrow(model_input_df))) {
  if (is.na(model_input_df$out_date[i])) {
    # right-hand censorship takes place
    cens_status[i]       <- 0
    outcome_cox_dates[i] <- model_input_df$end_date_outcome[i]
  } else {
    # event takes place (failure)
    cens_status[i]       <- 1
    outcome_cox_dates[i] <- model_input_df$out_date[i]
  }
}

# add data to dataframes
model_input_df$outcome_cox_dates <- as.numeric(outcome_cox_dates)
model_input_df$cens_status       <- cens_status

# calculate Nelson-Aalen estimator
H0              <- (survfit(Surv(outcome_cox_dates, cens_status) ~ 1, data = model_input_df) %>% summary(times = unique(model_input_df$outcome_cox_dates)))
H0              <- H0[c("time", "surv")]
names(H0)       <- c("outcome_cox_dates", "surv")
H0              <- as.data.frame(H0)
model_input_df <- merge(model_input_df, H0, all.x = TRUE, by = "outcome_cox_dates")
model_input_df <- rename(model_input_df, H0 = surv)


# Applying multiple imputation to BMI and smoking covariates for outcome -----
print("Applying multiple imputation to BMI and smoking covariates for outcome")

imp <- mice::mice(
  data       = model_input_df,
  m          = get_number_of_imputed_datasets(),
  maxit      = 20,
  formulas   = my_formulas,
  imp_method = unname(imp_method)
)

# Applying lasso, lasso_X and lasso_union models to the MI datasets
print("Applying lasso, lasso_X and lasso_union models to the MI datasets")

for (i in c(1:get_number_of_imputed_datasets())) {

  imp_data_i <- complete(imp, action = i)

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
    x      = lasso_cox_conf_matrix_preserving_factors_i,
    y      = lasso_cox_outcome_survival_i,
    nfolds = get_number_of_imputed_datasets(),
    family = "cox",      # cox regression
    alpha  = 1           # LASSO penalty
  )

  # tune regularisation parameter lambda to minimise cross-validated error (cvm)
  lambda_i <- cv_lasso_cox_model_i$lambda.min

  lasso_cox_model_i <- glmnet(
    x      = lasso_cox_conf_matrix_preserving_factors_i,
    y      = lasso_cox_outcome_survival_i,
    family = "cox",      # cox regression
    alpha  = 1,          # LASSO penalty
    lambda = lambda_i    # optimal lambda
  )

  lasso_cox_coefs_i        <- as.vector(lasso_cox_model_i$beta)
  names(lasso_cox_coefs_i) <- rownames(lasso_cox_model_i$beta)

  lasso_non_zero_vars_i  <- names(lasso_cox_coefs_i[lasso_cox_coefs_i != 0.0])
  lasso_vars_selected_i  <- convert_terms_to_vars(lasso_non_zero_vars_i)

  message("\n\ni")
  print(lasso_non_zero_vars_i)
}

stop("?")
