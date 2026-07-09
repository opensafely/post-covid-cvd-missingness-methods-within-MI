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


# Define unconfoundedness_test output folder -----------------------------------------
print("Creating output/unconfoundedness_test output folder")

unconfoundedness_test_dir <- "output/unconfoundedness_test/"
fs::dir_create(here::here(unconfoundedness_test_dir))


# Load data -------------------------------------------------------------------
print("Load data")

pooled_fully_adjusted_cox_model <- read.csv(
  paste0("output/all_cox_models/pooled_fully_adjusted_cox_model-", name, ".csv")
)

pooled_lasso_cox_model <- read.csv(
  paste0("output/all_cox_models/pooled_lasso_cox_model-", name, ".csv")
)

pooled_lasso_X_cox_model <- read.csv(
  paste0("output/all_cox_models/pooled_lasso_X_cox_model-", name, ".csv")
)

pooled_lasso_union_cox_model <- read.csv(
  paste0("output/all_cox_models/pooled_lasso_union_cox_model-", name, ".csv")
)

print(pooled_fully_adjusted_cox_model)
print(pooled_lasso_cox_model)
print(pooled_lasso_X_cox_model)
print(pooled_lasso_union_cox_model)
stop("?")
