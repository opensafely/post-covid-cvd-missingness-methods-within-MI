# ------------------------------------------------------------------------------
#
# unconfoundedness_test.R
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


# Define unconfoundedness_test output folder ----------------------------------
print("Creating output/unconfoundedness_test output folder")

unconfoundedness_test_dir <- "output/unconfoundedness_test/"
fs::dir_create(here::here(unconfoundedness_test_dir))



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

fully_adjusted_model_outcome_formula <- make_outcome_formula(vars_selected = fully_adjusted_var_selection_results)
lasso_model_outcome_formula          <- make_outcome_formula(vars_selected = lasso_var_selection_results)
lasso_X_model_outcome_formula        <- make_outcome_formula(vars_selected = lasso_X_var_selection_results)
lasso_union_model_outcome_formula    <- make_outcome_formula(vars_selected = lasso_union_var_selection_results)

fully_adjusted_model_exposure_formula <- make_exposure_formula(vars_selected = fully_adjusted_var_selection_results)
lasso_model_exposure_formula          <- make_exposure_formula(vars_selected = lasso_var_selection_results)
lasso_X_model_exposure_formula        <- make_exposure_formula(vars_selected = lasso_X_var_selection_results)
lasso_union_model_exposure_formula    <- make_exposure_formula(vars_selected = lasso_union_var_selection_results)


# All outcome regressions -------------------------------------------------
print("All outcome regressions")

# fully_adjusted
fully_adjusted_outcome_models <- with(
  data = imp,
  exp  = coxph(formula = as.formula(fully_adjusted_model_outcome_formula))
)
pooled_fully_adjusted_outcome_model <- summary(pool(fully_adjusted_outcome_models))

# lasso
lasso_outcome_models <- with(
  data = imp,
  exp  = coxph(formula = as.formula(lasso_model_outcome_formula))
)
pooled_lasso_outcome_model <- summary(pool(lasso_outcome_models))

# lasso_X
lasso_X_outcome_models <- with(
  data = imp,
  exp  = coxph(formula = as.formula(lasso_X_model_outcome_formula))
)
pooled_lasso_X_outcome_model <- summary(pool(lasso_X_outcome_models))

# lasso_union
lasso_union_outcome_models <- with(
  data = imp,
  exp  = coxph(formula = as.formula(lasso_union_model_outcome_formula))
)
pooled_lasso_union_outcome_model <- summary(pool(lasso_union_outcome_models))


# All exposure regressions -----------------------------------------------
print("All exposure regressions")

# fully_adjusted
fully_adjusted_exposure_models <- with(
  data = imp,
  exp  = glm(formula = as.formula(fully_adjusted_model_exposure_formula), family = "binomial")
)
pooled_fully_adjusted_exposure_model <- summary(pool(fully_adjusted_exposure_models))

# lasso
lasso_exposure_models <- with(
  data = imp,
  exp  = glm(formula = as.formula(lasso_model_exposure_formula), family = "binomial")
)
pooled_lasso_exposure_model <- summary(pool(lasso_exposure_models))

# lasso_X
lasso_X_exposure_models <- with(
  data = imp,
  exp  = glm(formula = as.formula(lasso_X_model_exposure_formula), family = "binomial")
)
pooled_lasso_X_exposure_model <- summary(pool(lasso_X_exposure_models))

# lasso_union
lasso_union_exposure_models <- with(
  data = imp,
  exp  = glm(formula = as.formula(lasso_union_model_exposure_formula), family = "binomial")
)
pooled_lasso_union_exposure_model <- summary(pool(lasso_union_exposure_models))


# Determine significant variables from regressions -----------------------
print("Determine significant variables from regressions")

fully_adjusted_outcome_model_significant_vars <- (
  pooled_fully_adjusted_outcome_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()

lasso_outcome_model_significant_vars <- (
  pooled_lasso_outcome_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()

lasso_X_outcome_model_significant_vars <- (
  pooled_lasso_X_outcome_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()

lasso_union_outcome_model_significant_vars <- (
  pooled_lasso_union_outcome_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()

fully_adjusted_exposure_model_significant_vars <- (
  pooled_fully_adjusted_exposure_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()

lasso_exposure_model_significant_vars <- (
  pooled_lasso_exposure_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()

lasso_X_exposure_model_significant_vars <- (
  pooled_lasso_X_exposure_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()

lasso_union_exposure_model_significant_vars <- (
  pooled_lasso_union_exposure_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()


pooled_fully_adjusted_outcome_model_significant_vars <- (
  pooled_fully_adjusted_outcome_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()

pooled_lasso_outcome_model_significant_vars <- (
  pooled_lasso_outcome_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()

pooled_lasso_X_outcome_model_significant_vars <- (
  pooled_lasso_X_outcome_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()

pooled_lasso_union_outcome_model_significant_vars <- (
  pooled_lasso_union_outcome_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()

pooled_fully_adjusted_exposure_model_significant_vars <- (
  pooled_fully_adjusted_exposure_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()

pooled_lasso_exposure_model_significant_vars <- (
  pooled_lasso_exposure_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()

pooled_lasso_X_exposure_model_significant_vars <- (
  pooled_lasso_X_exposure_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()

pooled_lasso_union_exposure_model_significant_vars <- (
  pooled_lasso_union_exposure_model %>% filter(p.value < 0.05)
)$term %>% as.vector() %>% convert_terms_to_vars()


# Determine test conditions for fully_adjusted set using above significance testing ---
print("Determine test conditions for fully_adjusted set using above significance testing")

# fully_adjusted_condition (i)
# Z is associated with (i.e., not independent of) X given all other covariates
fully_adjusted_condition_i        <- rep(FALSE, length.out = length(fully_adjusted_var_selection_results))
names(fully_adjusted_condition_i) <- fully_adjusted_var_selection_results
for (var in fully_adjusted_var_selection_results) {
  if (var %in% fully_adjusted_exposure_model_significant_vars) {
    fully_adjusted_condition_i[var] <- TRUE
  }
}

# fully_adjusted_condition (ii)
# Z and Y are fully_adjusted_conditionally independent given X and all other covariates
fully_adjusted_condition_ii        <- rep(FALSE, length.out = length(fully_adjusted_var_selection_results))
names(fully_adjusted_condition_ii) <- fully_adjusted_var_selection_results
for (var in fully_adjusted_var_selection_results) {
  if (!var %in% fully_adjusted_outcome_model_significant_vars) {
    fully_adjusted_condition_ii[var] <- TRUE
  }
}

# fully_adjusted_conditions (i) and (ii)
fully_adjusted_conditions_i_and_ii        <- fully_adjusted_condition_i & fully_adjusted_condition_ii
names(fully_adjusted_conditions_i_and_ii) <- fully_adjusted_var_selection_results

# test is TRUE if any covariate Z satisfies (i) and (ii)
# test is FALSE otherwise
fully_adjusted_conclusion <- any(fully_adjusted_conditions_i_and_ii)
fully_adjusted_conclusion_string <- ""
if (fully_adjusted_conclusion) {
  fully_adjusted_conclusion_string <- "Covariate set is sufficient for confounding adjustment"
} else {
  fully_adjusted_conclusion_string <- "Test is inconclusive, covariate set may or may not be sufficient"
}

fully_adjusted_test_table <- cbind(
  fully_adjusted_condition_i,
  fully_adjusted_condition_ii,
  fully_adjusted_conditions_i_and_ii
)
colnames(fully_adjusted_test_table) <- c("condition_i", "condition_ii", "condition_i_and_ii")
rownames(fully_adjusted_test_table) <- fully_adjusted_var_selection_results


# Determine test conditions for lasso set using above significance testing ---
print("Determine test conditions for lasso set using above significance testing")

# lasso_condition (i)
# Z is associated with (i.e., not independent of) X given all other covariates
lasso_condition_i        <- rep(FALSE, length.out = length(lasso_var_selection_results))
names(lasso_condition_i) <- lasso_var_selection_results
for (var in lasso_var_selection_results) {
  if (var %in% lasso_exposure_model_significant_vars) {
    lasso_condition_i[var] <- TRUE
  }
}

# lasso_condition (ii)
# Z and Y are lasso_conditionally independent given X and all other covariates
lasso_condition_ii        <- rep(FALSE, length.out = length(lasso_var_selection_results))
names(lasso_condition_ii) <- lasso_var_selection_results
for (var in lasso_var_selection_results) {
  if (!var %in% lasso_outcome_model_significant_vars) {
    lasso_condition_ii[var] <- TRUE
  }
}

# lasso_conditions (i) and (ii)
lasso_conditions_i_and_ii        <- lasso_condition_i & lasso_condition_ii
names(lasso_conditions_i_and_ii) <- lasso_var_selection_results

# test is TRUE if any covariate Z satisfies (i) and (ii)
# test is FALSE otherwise
lasso_conclusion <- any(lasso_conditions_i_and_ii)
lasso_conclusion_string <- ""
if (lasso_conclusion) {
  lasso_conclusion_string <- "Covariate set is sufficient for confounding adjustment"
} else {
  lasso_conclusion_string <- "Test is inconclusive, covariate set may or may not be sufficient"
}

lasso_test_table <- cbind(
  lasso_condition_i,
  lasso_condition_ii,
  lasso_conditions_i_and_ii
)
colnames(lasso_test_table) <- c("condition_i", "condition_ii", "condition_i_and_ii")
rownames(lasso_test_table) <- lasso_var_selection_results


# Determine test conditions for lasso_X set using above significance testing ---
print("Determine test conditions for lasso_X set using above significance testing")

# lasso_X_condition (i)
# Z is associated with (i.e., not independent of) X given all other covariates
lasso_X_condition_i        <- rep(FALSE, length.out = length(lasso_X_var_selection_results))
names(lasso_X_condition_i) <- lasso_X_var_selection_results
for (var in lasso_X_var_selection_results) {
  if (var %in% lasso_X_exposure_model_significant_vars) {
    lasso_X_condition_i[var] <- TRUE
  }
}

# lasso_X_condition (ii)
# Z and Y are lasso_X_conditionally independent given X and all other covariates
lasso_X_condition_ii        <- rep(FALSE, length.out = length(lasso_X_var_selection_results))
names(lasso_X_condition_ii) <- lasso_X_var_selection_results
for (var in lasso_X_var_selection_results) {
  if (!var %in% lasso_X_outcome_model_significant_vars) {
    lasso_X_condition_ii[var] <- TRUE
  }
}

# lasso_X_conditions (i) and (ii)
lasso_X_conditions_i_and_ii        <- lasso_X_condition_i & lasso_X_condition_ii
names(lasso_X_conditions_i_and_ii) <- lasso_X_var_selection_results

# test is TRUE if any covariate Z satisfies (i) and (ii)
# test is FALSE otherwise
lasso_X_conclusion <- any(lasso_X_conditions_i_and_ii)
lasso_X_conclusion_string <- ""
if (lasso_X_conclusion) {
  lasso_X_conclusion_string <- "Covariate set is sufficient for confounding adjustment"
} else {
  lasso_X_conclusion_string <- "Test is inconclusive, covariate set may or may not be sufficient"
}

lasso_X_test_table <- cbind(
  lasso_X_condition_i,
  lasso_X_condition_ii,
  lasso_X_conditions_i_and_ii
)
colnames(lasso_X_test_table) <- c("condition_i", "condition_ii", "condition_i_and_ii")
rownames(lasso_X_test_table) <- lasso_X_var_selection_results


# Determine test conditions for lasso_union set using above significance testing ---
print("Determine test conditions for lasso_union set using above significance testing")

# lasso_union_condition (i)
# Z is associated with (i.e., not independent of) X given all other covariates
lasso_union_condition_i        <- rep(FALSE, length.out = length(lasso_union_var_selection_results))
names(lasso_union_condition_i) <- lasso_union_var_selection_results
for (var in lasso_union_var_selection_results) {
  if (var %in% lasso_union_exposure_model_significant_vars) {
    lasso_union_condition_i[var] <- TRUE
  }
}

# lasso_union_condition (ii)
# Z and Y are lasso_union_conditionally independent given X and all other covariates
lasso_union_condition_ii        <- rep(FALSE, length.out = length(lasso_union_var_selection_results))
names(lasso_union_condition_ii) <- lasso_union_var_selection_results
for (var in lasso_union_var_selection_results) {
  if (!var %in% lasso_union_outcome_model_significant_vars) {
    lasso_union_condition_ii[var] <- TRUE
  }
}

# lasso_union_conditions (i) and (ii)
lasso_union_conditions_i_and_ii        <- lasso_union_condition_i & lasso_union_condition_ii
names(lasso_union_conditions_i_and_ii) <- lasso_union_var_selection_results

# test is TRUE if any covariate Z satisfies (i) and (ii)
# test is FALSE otherwise
lasso_union_conclusion <- any(lasso_union_conditions_i_and_ii)
lasso_union_conclusion_string <- ""
if (lasso_union_conclusion) {
  lasso_union_conclusion_string <- "Covariate set is sufficient for confounding adjustment"
} else {
  lasso_union_conclusion_string <- "Test is inconclusive, covariate set may or may not be sufficient"
}

lasso_union_test_table <- cbind(
  lasso_union_condition_i,
  lasso_union_condition_ii,
  lasso_union_conditions_i_and_ii
)
colnames(lasso_union_test_table) <- c("condition_i", "condition_ii", "condition_i_and_ii")
rownames(lasso_union_test_table) <- lasso_union_var_selection_results


# Conclusion table for all var sets -------------------------------------------
print("Conclusion table for all var sets")

all_var_sets_conclusion_table <- cbind(
  c("fully_adjusted", "lasso", "lasso_X", "lasso_union"),
  c(fully_adjusted_conclusion, lasso_conclusion, lasso_X_conclusion, lasso_union_conclusion),
  c(fully_adjusted_conclusion_string, lasso_conclusion_string, lasso_X_conclusion_string, lasso_union_conclusion_string)
)


# Save results ----------------------------------------------------------------

write.csv(
  all_var_sets_conclusion_table,
  paste0(unconfoundedness_test_dir, "all_var_sets_conclusion_table-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  pooled_fully_adjusted_outcome_model,
  paste0(unconfoundedness_test_dir, "fully_adjusted_outcome_model_results-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  pooled_fully_adjusted_exposure_model,
  paste0(unconfoundedness_test_dir, "fully_adjusted_exposure_model_results-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  fully_adjusted_test_table,
  paste0(unconfoundedness_test_dir, "fully_adjusted_test_table-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  pooled_lasso_outcome_model,
  paste0(unconfoundedness_test_dir, "lasso_outcome_model_results-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  pooled_lasso_exposure_model,
  paste0(unconfoundedness_test_dir, "lasso_exposure_model_results-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  lasso_test_table,
  paste0(unconfoundedness_test_dir, "lasso_test_table-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  pooled_lasso_X_outcome_model,
  paste0(unconfoundedness_test_dir, "lasso_X_outcome_model_results-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  pooled_lasso_X_exposure_model,
  paste0(unconfoundedness_test_dir, "lasso_X_exposure_model_results-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  lasso_X_test_table,
  paste0(unconfoundedness_test_dir, "lasso_X_test_table-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  pooled_lasso_union_outcome_model,
  paste0(unconfoundedness_test_dir, "lasso_union_outcome_model_results-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  pooled_lasso_union_exposure_model,
  paste0(unconfoundedness_test_dir, "lasso_union_exposure_model_results-", name, ".csv"),
  row.names = TRUE
)

write.csv(
  lasso_union_test_table,
  paste0(unconfoundedness_test_dir, "lasso_union_test_table-", name, ".csv"),
  row.names = TRUE
)
