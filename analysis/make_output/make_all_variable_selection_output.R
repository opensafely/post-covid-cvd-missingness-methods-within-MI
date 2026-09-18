# ------------------------------------------------------------------------------
#
# make_all_variable_selection_output.R
#
# Tidy all outputs from the all lasso, lasso_X and lasso_union models
#
# Arguments:
#  - None!
#
# Returns:
#  - All lasso, lasso_X and lasso_union model outputs
#
# Authors: Emma Tarmey
#
# ------------------------------------------------------------------------------


# Load packages ----------------------------------------------------------------
print('Load packages')

library(magrittr)
library(data.table)
library(stringr)
library(tidyr)


# Source common functions ------------------------------------------------------
print('Source common functions')

source("analysis/utility.R")


# Define make_output folder ----------------------------------------------------
print("Creating output/make_output output folder")

makeout_dir <- "output/make_output/"
fs::dir_create(here::here(makeout_dir))


# Load all cox models files -----------------------------------------------
print("Load all cox models files")


lasso_cohort_prevax_main_ami                             <- read.csv("output/all_variable_selection/lasso_aggregate_var_selection_results-cohort_prevax-main-ami.csv")
lasso_cohort_prevax_main_stroke_sahhs                    <- read.csv("output/all_variable_selection/lasso_aggregate_var_selection_results-cohort_prevax-main-stroke_sahhs.csv")
lasso_cohort_prevax_sub_covidhospital_FALSE_ami          <- read.csv("output/all_variable_selection/lasso_aggregate_var_selection_results-cohort_prevax-sub_covidhospital_FALSE-ami.csv")
lasso_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- read.csv("output/all_variable_selection/lasso_aggregate_var_selection_results-cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs.csv")
lasso_cohort_prevax_sub_covidhospital_TRUE_ami           <- read.csv("output/all_variable_selection/lasso_aggregate_var_selection_results-cohort_prevax-sub_covidhospital_TRUE-ami.csv")
lasso_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs  <- read.csv("output/all_variable_selection/lasso_aggregate_var_selection_results-cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs.csv")

lasso_X_cohort_prevax_main_ami                             <- read.csv("output/all_variable_selection/lasso_X_aggregate_var_selection_results-cohort_prevax-main-ami.csv")
lasso_X_cohort_prevax_main_stroke_sahhs                    <- read.csv("output/all_variable_selection/lasso_X_aggregate_var_selection_results-cohort_prevax-main-stroke_sahhs.csv")
lasso_X_cohort_prevax_sub_covidhospital_FALSE_ami          <- read.csv("output/all_variable_selection/lasso_X_aggregate_var_selection_results-cohort_prevax-sub_covidhospital_FALSE-ami.csv")
lasso_X_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- read.csv("output/all_variable_selection/lasso_X_aggregate_var_selection_results-cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs.csv")
lasso_X_cohort_prevax_sub_covidhospital_TRUE_ami           <- read.csv("output/all_variable_selection/lasso_X_aggregate_var_selection_results-cohort_prevax-sub_covidhospital_TRUE-ami.csv")
lasso_X_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs  <- read.csv("output/all_variable_selection/lasso_X_aggregate_var_selection_results-cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs.csv")

lasso_union_cohort_prevax_main_ami                             <- read.csv("output/all_variable_selection/lasso_union_aggregate_var_selection_results-cohort_prevax-main-ami.csv")
lasso_union_cohort_prevax_main_stroke_sahhs                    <- read.csv("output/all_variable_selection/lasso_union_aggregate_var_selection_results-cohort_prevax-main-stroke_sahhs.csv")
lasso_union_cohort_prevax_sub_covidhospital_FALSE_ami          <- read.csv("output/all_variable_selection/lasso_union_aggregate_var_selection_results-cohort_prevax-sub_covidhospital_FALSE-ami.csv")
lasso_union_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- read.csv("output/all_variable_selection/lasso_union_aggregate_var_selection_results-cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs.csv")
lasso_union_cohort_prevax_sub_covidhospital_TRUE_ami           <- read.csv("output/all_variable_selection/lasso_union_aggregate_var_selection_results-cohort_prevax-sub_covidhospital_TRUE-ami.csv")
lasso_union_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs  <- read.csv("output/all_variable_selection/lasso_union_aggregate_var_selection_results-cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs.csv")

lasso_mean_cohort_prevax_main_ami                             <- read.csv("output/all_variable_selection/lasso_mean_var_selection_results-cohort_prevax-main-ami.csv")
lasso_mean_cohort_prevax_main_stroke_sahhs                    <- read.csv("output/all_variable_selection/lasso_mean_var_selection_results-cohort_prevax-main-stroke_sahhs.csv")
lasso_mean_cohort_prevax_sub_covidhospital_FALSE_ami          <- read.csv("output/all_variable_selection/lasso_mean_var_selection_results-cohort_prevax-sub_covidhospital_FALSE-ami.csv")
lasso_mean_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- read.csv("output/all_variable_selection/lasso_mean_var_selection_results-cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs.csv")
lasso_mean_cohort_prevax_sub_covidhospital_TRUE_ami           <- read.csv("output/all_variable_selection/lasso_mean_var_selection_results-cohort_prevax-sub_covidhospital_TRUE-ami.csv")
lasso_mean_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs  <- read.csv("output/all_variable_selection/lasso_mean_var_selection_results-cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs.csv")

lasso_X_mean_cohort_prevax_main_ami                             <- read.csv("output/all_variable_selection/lasso_X_mean_var_selection_results-cohort_prevax-main-ami.csv")
lasso_X_mean_cohort_prevax_main_stroke_sahhs                    <- read.csv("output/all_variable_selection/lasso_X_mean_var_selection_results-cohort_prevax-main-stroke_sahhs.csv")
lasso_X_mean_cohort_prevax_sub_covidhospital_FALSE_ami          <- read.csv("output/all_variable_selection/lasso_X_mean_var_selection_results-cohort_prevax-sub_covidhospital_FALSE-ami.csv")
lasso_X_mean_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- read.csv("output/all_variable_selection/lasso_X_mean_var_selection_results-cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs.csv")
lasso_X_mean_cohort_prevax_sub_covidhospital_TRUE_ami           <- read.csv("output/all_variable_selection/lasso_X_mean_var_selection_results-cohort_prevax-sub_covidhospital_TRUE-ami.csv")
lasso_X_mean_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs  <- read.csv("output/all_variable_selection/lasso_X_mean_var_selection_results-cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs.csv")

lasso_union_mean_cohort_prevax_main_ami                             <- read.csv("output/all_variable_selection/lasso_union_mean_var_selection_results-cohort_prevax-main-ami.csv")
lasso_union_mean_cohort_prevax_main_stroke_sahhs                    <- read.csv("output/all_variable_selection/lasso_union_mean_var_selection_results-cohort_prevax-main-stroke_sahhs.csv")
lasso_union_mean_cohort_prevax_sub_covidhospital_FALSE_ami          <- read.csv("output/all_variable_selection/lasso_union_mean_var_selection_results-cohort_prevax-sub_covidhospital_FALSE-ami.csv")
lasso_union_mean_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- read.csv("output/all_variable_selection/lasso_union_mean_var_selection_results-cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs.csv")
lasso_union_mean_cohort_prevax_sub_covidhospital_TRUE_ami           <- read.csv("output/all_variable_selection/lasso_union_mean_var_selection_results-cohort_prevax-sub_covidhospital_TRUE-ami.csv")
lasso_union_mean_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs  <- read.csv("output/all_variable_selection/lasso_union_mean_var_selection_results-cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs.csv")


# Add method and name columns to each lasso table ----------------------
print("Add method and name columns to each lasso table")

lasso_cohort_prevax_main_ami <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-main-ami",
  lasso_cohort_prevax_main_ami
)

lasso_cohort_prevax_main_stroke_sahhs <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-main-stroke_sahhs",
  lasso_cohort_prevax_main_stroke_sahhs
)

lasso_cohort_prevax_sub_covidhospital_FALSE_ami <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-sub_covidhospital_FALSE-ami",
  lasso_cohort_prevax_sub_covidhospital_FALSE_ami
)

lasso_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs",
  lasso_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs
)

lasso_cohort_prevax_sub_covidhospital_TRUE_ami <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-sub_covidhospital_TRUE-ami",
  lasso_cohort_prevax_sub_covidhospital_TRUE_ami
)

lasso_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs",
  lasso_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Stack lasso tables -----
print("Stack lasso tables")

lasso_all_aggregate_outputs <- rbind(
  lasso_cohort_prevax_main_ami,
  lasso_cohort_prevax_main_stroke_sahhs,
  lasso_cohort_prevax_sub_covidhospital_FALSE_ami,
  lasso_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs,
  lasso_cohort_prevax_sub_covidhospital_TRUE_ami,
  lasso_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Add method and name columns to each lasso_X table ----------------------
print("Add method and name columns to each lasso_X table")

lasso_X_cohort_prevax_main_ami <- cbind(
  method   = "lasso_X",
  name     = "cohort_prevax-main-ami",
  lasso_X_cohort_prevax_main_ami
)

lasso_X_cohort_prevax_main_stroke_sahhs <- cbind(
  method   = "lasso_X",
  name     = "cohort_prevax-main-stroke_sahhs",
  lasso_X_cohort_prevax_main_stroke_sahhs
)

lasso_X_cohort_prevax_sub_covidhospital_FALSE_ami <- cbind(
  method   = "lasso_X",
  name     = "cohort_prevax-sub_covidhospital_FALSE-ami",
  lasso_X_cohort_prevax_sub_covidhospital_FALSE_ami
)

lasso_X_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- cbind(
  method   = "lasso_X",
  name     = "cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs",
  lasso_X_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs
)

lasso_X_cohort_prevax_sub_covidhospital_TRUE_ami <- cbind(
  method   = "lasso_X",
  name     = "cohort_prevax-sub_covidhospital_TRUE-ami",
  lasso_X_cohort_prevax_sub_covidhospital_TRUE_ami
)

lasso_X_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs <- cbind(
  method   = "lasso_X",
  name     = "cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs",
  lasso_X_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Stack lasso_X tables -----
print("Stack lasso_X tables")

lasso_X_all_aggregate_outputs <- rbind(
  lasso_X_cohort_prevax_main_ami,
  lasso_X_cohort_prevax_main_stroke_sahhs,
  lasso_X_cohort_prevax_sub_covidhospital_FALSE_ami,
  lasso_X_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs,
  lasso_X_cohort_prevax_sub_covidhospital_TRUE_ami,
  lasso_X_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Add method and name columns to each lasso_union table ----------------------
print("Add method and name columns to each lasso_union table")

lasso_union_cohort_prevax_main_ami <- cbind(
  method   = "lasso_union",
  name     = "cohort_prevax-main-ami",
  lasso_union_cohort_prevax_main_ami
)

lasso_union_cohort_prevax_main_stroke_sahhs <- cbind(
  method   = "lasso_union",
  name     = "cohort_prevax-main-stroke_sahhs",
  lasso_union_cohort_prevax_main_stroke_sahhs
)

lasso_union_cohort_prevax_sub_covidhospital_FALSE_ami <- cbind(
  method   = "lasso_union",
  name     = "cohort_prevax-sub_covidhospital_FALSE-ami",
  lasso_union_cohort_prevax_sub_covidhospital_FALSE_ami
)

lasso_union_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- cbind(
  method   = "lasso_union",
  name     = "cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs",
  lasso_union_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs
)

lasso_union_cohort_prevax_sub_covidhospital_TRUE_ami <- cbind(
  method   = "lasso_union",
  name     = "cohort_prevax-sub_covidhospital_TRUE-ami",
  lasso_union_cohort_prevax_sub_covidhospital_TRUE_ami
)

lasso_union_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs <- cbind(
  method   = "lasso_union",
  name     = "cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs",
  lasso_union_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Stack lasso_union tables -----
print("Stack lasso_union tables")

lasso_union_all_aggregate_outputs <- rbind(
  lasso_union_cohort_prevax_main_ami,
  lasso_union_cohort_prevax_main_stroke_sahhs,
  lasso_union_cohort_prevax_sub_covidhospital_FALSE_ami,
  lasso_union_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs,
  lasso_union_cohort_prevax_sub_covidhospital_TRUE_ami,
  lasso_union_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Add method and name columns to each lasso table ----------------------
print("Add method and name columns to each lasso table")

lasso_mean_cohort_prevax_main_ami <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-main-ami",
  lasso_mean_cohort_prevax_main_ami
)

lasso_mean_cohort_prevax_main_stroke_sahhs <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-main-stroke_sahhs",
  lasso_mean_cohort_prevax_main_stroke_sahhs
)

lasso_mean_cohort_prevax_sub_covidhospital_FALSE_ami <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-sub_covidhospital_FALSE-ami",
  lasso_mean_cohort_prevax_sub_covidhospital_FALSE_ami
)

lasso_mean_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs",
  lasso_mean_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs
)

lasso_mean_cohort_prevax_sub_covidhospital_TRUE_ami <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-sub_covidhospital_TRUE-ami",
  lasso_mean_cohort_prevax_sub_covidhospital_TRUE_ami
)

lasso_mean_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs",
  lasso_mean_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Stack lasso tables -----
print("Stack lasso tables")

lasso_all_mean_outputs <- rbind(
  lasso_mean_cohort_prevax_main_ami,
  lasso_mean_cohort_prevax_main_stroke_sahhs,
  lasso_mean_cohort_prevax_sub_covidhospital_FALSE_ami,
  lasso_mean_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs,
  lasso_mean_cohort_prevax_sub_covidhospital_TRUE_ami,
  lasso_mean_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Add method and name columns to each lasso_X_mean table ----------------------
print("Add method and name columns to each lasso_X_mean table")

lasso_X_mean_cohort_prevax_main_ami <- cbind(
  method   = "lasso_X_mean",
  name     = "cohort_prevax-main-ami",
  lasso_X_mean_cohort_prevax_main_ami
)

lasso_X_mean_cohort_prevax_main_stroke_sahhs <- cbind(
  method   = "lasso_X_mean",
  name     = "cohort_prevax-main-stroke_sahhs",
  lasso_X_mean_cohort_prevax_main_stroke_sahhs
)

lasso_X_mean_cohort_prevax_sub_covidhospital_FALSE_ami <- cbind(
  method   = "lasso_X_mean",
  name     = "cohort_prevax-sub_covidhospital_FALSE-ami",
  lasso_X_mean_cohort_prevax_sub_covidhospital_FALSE_ami
)

lasso_X_mean_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- cbind(
  method   = "lasso_X_mean",
  name     = "cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs",
  lasso_X_mean_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs
)

lasso_X_mean_cohort_prevax_sub_covidhospital_TRUE_ami <- cbind(
  method   = "lasso_X_mean",
  name     = "cohort_prevax-sub_covidhospital_TRUE-ami",
  lasso_X_mean_cohort_prevax_sub_covidhospital_TRUE_ami
)

lasso_X_mean_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs <- cbind(
  method   = "lasso_X_mean",
  name     = "cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs",
  lasso_X_mean_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Stack lasso_X_mean tables -----
print("Stack lasso_X_mean tables")

lasso_X_all_mean_outputs <- rbind(
  lasso_X_mean_cohort_prevax_main_ami,
  lasso_X_mean_cohort_prevax_main_stroke_sahhs,
  lasso_X_mean_cohort_prevax_sub_covidhospital_FALSE_ami,
  lasso_X_mean_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs,
  lasso_X_mean_cohort_prevax_sub_covidhospital_TRUE_ami,
  lasso_X_mean_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Add method and name columns to each lasso_union_mean table ----------------------
print("Add method and name columns to each lasso_union_mean table")

lasso_union_mean_cohort_prevax_main_ami <- cbind(
  method   = "lasso_union_mean",
  name     = "cohort_prevax-main-ami",
  lasso_union_mean_cohort_prevax_main_ami
)

lasso_union_mean_cohort_prevax_main_stroke_sahhs <- cbind(
  method   = "lasso_union_mean",
  name     = "cohort_prevax-main-stroke_sahhs",
  lasso_union_mean_cohort_prevax_main_stroke_sahhs
)

lasso_union_mean_cohort_prevax_sub_covidhospital_FALSE_ami <- cbind(
  method   = "lasso_union_mean",
  name     = "cohort_prevax-sub_covidhospital_FALSE-ami",
  lasso_union_mean_cohort_prevax_sub_covidhospital_FALSE_ami
)

lasso_union_mean_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- cbind(
  method   = "lasso_union_mean",
  name     = "cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs",
  lasso_union_mean_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs
)

lasso_union_mean_cohort_prevax_sub_covidhospital_TRUE_ami <- cbind(
  method   = "lasso_union_mean",
  name     = "cohort_prevax-sub_covidhospital_TRUE-ami",
  lasso_union_mean_cohort_prevax_sub_covidhospital_TRUE_ami
)

lasso_union_mean_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs <- cbind(
  method   = "lasso_union_mean",
  name     = "cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs",
  lasso_union_mean_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Stack lasso_union_mean tables -----
print("Stack lasso_union_mean tables")

lasso_union_all_mean_outputs <- rbind(
  lasso_union_mean_cohort_prevax_main_ami,
  lasso_union_mean_cohort_prevax_main_stroke_sahhs,
  lasso_union_mean_cohort_prevax_sub_covidhospital_FALSE_ami,
  lasso_union_mean_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs,
  lasso_union_mean_cohort_prevax_sub_covidhospital_TRUE_ami,
  lasso_union_mean_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Save all tables -----
print("Save all tables")

all_aggregate_outputs <- rbind(
  lasso_all_aggregate_outputs,
  lasso_X_all_aggregate_outputs,
  lasso_union_all_aggregate_outputs
)
colnames(all_aggregate_outputs) <- c("method", "name", "covariates")

all_mean_outputs <- rbind(
  lasso_all_mean_outputs,
  lasso_X_all_mean_outputs,
  lasso_union_all_mean_outputs
)
colnames(all_mean_outputs) <- c("method", "name", "covariates", "mean")


write.csv(
  all_aggregate_outputs,
  paste0(makeout_dir, "all_aggregate_variable_selection_outputs.csv"),
  row.names = FALSE
)

write.csv(
  all_mean_outputs,
  paste0(makeout_dir, "all_mean_variable_selection_outputs.csv"),
  row.names = FALSE
)
