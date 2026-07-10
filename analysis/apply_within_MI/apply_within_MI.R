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

# drop now unused factor level "Missing"
model_input_df$cov_cat_smoking <- factor(
  model_input_df$cov_cat_smoking,
  levels = levels(droplevels(model_input_df$cov_cat_smoking))
)

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

all_var_names_except_bmi <- c(
  "cov_bin_ami",
  "cov_bin_sahhs",
  "cov_bin_covid",

  "cov_num_age",
  "cov_cat_sex",
  # "cov_num_bmi",
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

all_var_names_except_smoking <- c(
  "cov_bin_ami",
  "cov_bin_sahhs",
  "cov_bin_covid",

  "cov_num_age",
  "cov_cat_sex",
  "cov_num_bmi",
  "cov_cat_ethnicity",
  "cov_cat_imd",
  # "cov_cat_smoking",

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

my_formulas <- list(
  cov_cat_smoking = as.formula(paste0("cov_cat_smoking ~ ", paste(all_var_names_except_smoking, collapse = " + "), " + H0")),
  cov_num_bmi     = as.formula(paste0("cov_num_bmi ~ ",     paste(all_var_names_except_bmi, collapse = " + "), " + H0"))
)


# Calculate Nelson-Aalen Estimator for outcome -----------
print("Calculate Nelson-Aalen Estimator for outcome")

basehaz_with_SE <- function(fit, newdata, centered = TRUE) 
{
  if (inherits(fit, "coxphms"))
    stop("the basehaz function is not implemented for multi-state models")
  if (!inherits(fit, "coxph")) 
    stop("must be a coxph object")
  if (!missing(newdata)) {
    sfit <- survfit(fit, newdata=newdata, se.fit=TRUE)
    chaz <- sfit$cumhaz
  }
  else {
    sfit <- survfit(fit, se.fit=TRUE)
    if (!centered) {
      # The right thing to do here is to call survfit with a vector of
      #  all zeros for the "subject to predict".  But if there is a factor
      #  in the model, there may be no subject at all who will give all
      #  zeros, so we post process instead
      zcoef <- ifelse(is.na(coef(fit)), 0, coef(fit))
      offset <- sum(fit$means * zcoef)
      chaz <- sfit$cumhaz * exp(-offset)
    }
    else {
      chaz <- sfit$cumhaz
    }
  }

  new <- data.frame(
    hazard         = chaz,
    time           = sfit$time,
    std_err_cumhaz = sfit$std.err
  )

  strata <- sfit$strata
  if (!is.null(strata)) {
    new$strata <- factor(rep(names(strata), strata), levels = names(strata))
  }

  return (new)
}

nelsonaalen_with_SE <- function(data, timevar, statusvar, ...) {
  if (!is.data.frame(data)) {
    stop("Data must be a data frame")
  }
  timevar <- as.character(substitute(timevar))
  statusvar <- as.character(substitute(statusvar))
  time <- data[, timevar, drop = TRUE]
  status <- data[, statusvar, drop = TRUE]

  coxph_obj <- survival::coxph(survival::Surv(time, status) ~ 1, ...)
  hazard <- basehaz_with_SE(coxph_obj)

  # Adjust depending on near-tie correction
  idx <- if (coxph_obj$timefix) {
    match(coxph_obj$y[, "time"], hazard[, "time"])
  } else match(time, hazard[, "time"])
  
  nelsonaalen_estimates_with_SE <- data.frame(
    nelsonaalen_estimates = hazard[idx, "hazard"],
    nelsonaalen_se        = hazard[idx, "std_err_cumhaz"]
  )

  return (nelsonaalen_estimates_with_SE)
}

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
model_input_df_nelsonaalen    <- data.frame(time = model_input_df$outcome_cox_dates, status = model_input_df$cens_status)
H0                            <- nelsonaalen_with_SE(model_input_df_nelsonaalen, time, status)
model_input_df_nelsonaalen$H0 <- H0$nelsonaalen_estimates
model_input_df_nelsonaalen$se <- H0$nelsonaalen_se
model_input_df$H0             <- H0$nelsonaalen_estimates


# Applying multiple imputation to BMI and smoking covariates for outcome -----
print("Applying multiple imputation to BMI and smoking covariates for outcome")

imp <- mice::mice(
  data       = model_input_df,
  m          = get_number_of_imputed_datasets(),
  maxit      = 20,
  formulas   = my_formulas,
  imp_method = unname(imp_method)
)

saveRDS(
  imp,
  paste0(apply_within_MI_dir, "apply_within_MI_imp_datasets_", name, ".rds")
)
