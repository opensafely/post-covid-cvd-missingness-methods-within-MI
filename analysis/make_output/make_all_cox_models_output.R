# ------------------------------------------------------------------------------
#
# make_all_cox_models_output.R
#
# Tidy all outputs from the all cox models
#
# Arguments:
#  - None!
#
# Returns:
#  - All cix model outputs
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

fully_adjusted_cohort_prevax_main_ami                             <- read.csv("output/all_cox_models/fully_adjusted_pooled_cox_results-cohort_prevax-main-ami.csv")
fully_adjusted_cohort_prevax_main_stroke_sahhs                    <- read.csv("output/all_cox_models/fully_adjusted_pooled_cox_results-cohort_prevax-main-stroke_sahhs.csv")
fully_adjusted_cohort_prevax_sub_covidhospital_FALSE_ami          <- read.csv("output/all_cox_models/fully_adjusted_pooled_cox_results-cohort_prevax-sub_covidhospital_FALSE-ami.csv")
fully_adjusted_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- read.csv("output/all_cox_models/fully_adjusted_pooled_cox_results-cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs.csv")
fully_adjusted_cohort_prevax_sub_covidhospital_TRUE_ami           <- read.csv("output/all_cox_models/fully_adjusted_pooled_cox_results-cohort_prevax-sub_covidhospital_TRUE-ami.csv")
fully_adjusted_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs  <- read.csv("output/all_cox_models/fully_adjusted_pooled_cox_results-cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs.csv")

lasso_cohort_prevax_main_ami                             <- read.csv("output/all_cox_models/lasso_pooled_cox_results-cohort_prevax-main-ami.csv")
lasso_cohort_prevax_main_stroke_sahhs                    <- read.csv("output/all_cox_models/lasso_pooled_cox_results-cohort_prevax-main-stroke_sahhs.csv")
lasso_cohort_prevax_sub_covidhospital_FALSE_ami          <- read.csv("output/all_cox_models/lasso_pooled_cox_results-cohort_prevax-sub_covidhospital_FALSE-ami.csv")
lasso_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- read.csv("output/all_cox_models/lasso_pooled_cox_results-cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs.csv")
lasso_cohort_prevax_sub_covidhospital_TRUE_ami           <- read.csv("output/all_cox_models/lasso_pooled_cox_results-cohort_prevax-sub_covidhospital_TRUE-ami.csv")
lasso_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs  <- read.csv("output/all_cox_models/lasso_pooled_cox_results-cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs.csv")

lasso_X_cohort_prevax_main_ami                             <- read.csv("output/all_cox_models/lasso_X_pooled_cox_results-cohort_prevax-main-ami.csv")
lasso_X_cohort_prevax_main_stroke_sahhs                    <- read.csv("output/all_cox_models/lasso_X_pooled_cox_results-cohort_prevax-main-stroke_sahhs.csv")
lasso_X_cohort_prevax_sub_covidhospital_FALSE_ami          <- read.csv("output/all_cox_models/lasso_X_pooled_cox_results-cohort_prevax-sub_covidhospital_FALSE-ami.csv")
lasso_X_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- read.csv("output/all_cox_models/lasso_X_pooled_cox_results-cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs.csv")
lasso_X_cohort_prevax_sub_covidhospital_TRUE_ami           <- read.csv("output/all_cox_models/lasso_X_pooled_cox_results-cohort_prevax-sub_covidhospital_TRUE-ami.csv")
lasso_X_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs  <- read.csv("output/all_cox_models/lasso_X_pooled_cox_results-cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs.csv")

lasso_union_cohort_prevax_main_ami                             <- read.csv("output/all_cox_models/lasso_union_pooled_cox_results-cohort_prevax-main-ami.csv")
lasso_union_cohort_prevax_main_stroke_sahhs                    <- read.csv("output/all_cox_models/lasso_union_pooled_cox_results-cohort_prevax-main-stroke_sahhs.csv")
lasso_union_cohort_prevax_sub_covidhospital_FALSE_ami          <- read.csv("output/all_cox_models/lasso_union_pooled_cox_results-cohort_prevax-sub_covidhospital_FALSE-ami.csv")
lasso_union_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- read.csv("output/all_cox_models/lasso_union_pooled_cox_results-cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs.csv")
lasso_union_cohort_prevax_sub_covidhospital_TRUE_ami           <- read.csv("output/all_cox_models/lasso_union_pooled_cox_results-cohort_prevax-sub_covidhospital_TRUE-ami.csv")
lasso_union_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs  <- read.csv("output/all_cox_models/lasso_union_pooled_cox_results-cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs.csv")


# Add method and name columns to each fully_adjusted cox table ----------------------
print("Add method and name columns to each fully_adjusted cox table")

fully_adjusted_cohort_prevax_main_ami <- cbind(
  method   = "fully_adjusted",
  name     = "cohort_prevax-main-ami",
  fully_adjusted_cohort_prevax_main_ami
)

fully_adjusted_cohort_prevax_main_ami <- subset(
  fully_adjusted_cohort_prevax_main_ami,
  select = -c(X)
)

fully_adjusted_cohort_prevax_main_stroke_sahhs <- cbind(
  method   = "fully_adjusted",
  name     = "cohort_prevax-main-stroke_sahhs",
  fully_adjusted_cohort_prevax_main_stroke_sahhs
)

fully_adjusted_cohort_prevax_main_stroke_sahhs <- subset(
  fully_adjusted_cohort_prevax_main_stroke_sahhs,
  select = -c(X)
)

fully_adjusted_cohort_prevax_sub_covidhospital_FALSE_ami <- cbind(
  method   = "fully_adjusted",
  name     = "cohort_prevax-sub_covidhospital_FALSE-ami",
  fully_adjusted_cohort_prevax_sub_covidhospital_FALSE_ami
)

fully_adjusted_cohort_prevax_sub_covidhospital_FALSE_ami <- subset(
  fully_adjusted_cohort_prevax_sub_covidhospital_FALSE_ami,
  select = -c(X)
)

fully_adjusted_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- cbind(
  method   = "fully_adjusted",
  name     = "cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs",
  fully_adjusted_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs
)

fully_adjusted_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- subset(
  fully_adjusted_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs,
  select = -c(X)
)

fully_adjusted_cohort_prevax_sub_covidhospital_TRUE_ami <- cbind(
  method   = "fully_adjusted",
  name     = "cohort_prevax-sub_covidhospital_TRUE-ami",
  fully_adjusted_cohort_prevax_sub_covidhospital_TRUE_ami
)

fully_adjusted_cohort_prevax_sub_covidhospital_TRUE_ami <- subset(
  fully_adjusted_cohort_prevax_sub_covidhospital_TRUE_ami,
  select = -c(X)
)

fully_adjusted_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs <- cbind(
  method   = "fully_adjusted",
  name     = "cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs",
  fully_adjusted_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)

fully_adjusted_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs <- subset(
  fully_adjusted_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs,
  select = -c(X)
)


# Stack fully_adjusted tables -----
print("Stack fully_adjusted tables")

fully_adjusted_all_cox_outputs <- rbind(
  fully_adjusted_cohort_prevax_main_ami,
  fully_adjusted_cohort_prevax_main_stroke_sahhs,
  fully_adjusted_cohort_prevax_sub_covidhospital_FALSE_ami,
  fully_adjusted_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs,
  fully_adjusted_cohort_prevax_sub_covidhospital_TRUE_ami,
  fully_adjusted_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Add method and name columns to each lasso cox table ----------------------
print("Add method and name columns to each lasso cox table")

lasso_cohort_prevax_main_ami <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-main-ami",
  lasso_cohort_prevax_main_ami
)

lasso_cohort_prevax_main_ami <- subset(
  lasso_cohort_prevax_main_ami,
  select = -c(X)
)

lasso_cohort_prevax_main_stroke_sahhs <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-main-stroke_sahhs",
  lasso_cohort_prevax_main_stroke_sahhs
)

lasso_cohort_prevax_main_stroke_sahhs <- subset(
  lasso_cohort_prevax_main_stroke_sahhs,
  select = -c(X)
)

lasso_cohort_prevax_sub_covidhospital_FALSE_ami <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-sub_covidhospital_FALSE-ami",
  lasso_cohort_prevax_sub_covidhospital_FALSE_ami
)

lasso_cohort_prevax_sub_covidhospital_FALSE_ami <- subset(
  lasso_cohort_prevax_sub_covidhospital_FALSE_ami,
  select = -c(X)
)

lasso_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs",
  lasso_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs
)

lasso_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- subset(
  lasso_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs,
  select = -c(X)
)

lasso_cohort_prevax_sub_covidhospital_TRUE_ami <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-sub_covidhospital_TRUE-ami",
  lasso_cohort_prevax_sub_covidhospital_TRUE_ami
)

lasso_cohort_prevax_sub_covidhospital_TRUE_ami <- subset(
  lasso_cohort_prevax_sub_covidhospital_TRUE_ami,
  select = -c(X)
)

lasso_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs <- cbind(
  method   = "lasso",
  name     = "cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs",
  lasso_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)

lasso_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs <- subset(
  lasso_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs,
  select = -c(X)
)


# Stack lasso tables -----
print("Stack lasso tables")

lasso_all_cox_outputs <- rbind(
  lasso_cohort_prevax_main_ami,
  lasso_cohort_prevax_main_stroke_sahhs,
  lasso_cohort_prevax_sub_covidhospital_FALSE_ami,
  lasso_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs,
  lasso_cohort_prevax_sub_covidhospital_TRUE_ami,
  lasso_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Add method and name columns to each lasso_X cox table ----------------------
print("Add method and name columns to each lasso_X cox table")

lasso_X_cohort_prevax_main_ami <- cbind(
  method   = "lasso_X",
  name     = "cohort_prevax-main-ami",
  lasso_X_cohort_prevax_main_ami
)

lasso_X_cohort_prevax_main_ami <- subset(
  lasso_X_cohort_prevax_main_ami,
  select = -c(X)
)

lasso_X_cohort_prevax_main_stroke_sahhs <- cbind(
  method   = "lasso_X",
  name     = "cohort_prevax-main-stroke_sahhs",
  lasso_X_cohort_prevax_main_stroke_sahhs
)

lasso_X_cohort_prevax_main_stroke_sahhs <- subset(
  lasso_X_cohort_prevax_main_stroke_sahhs,
  select = -c(X)
)

lasso_X_cohort_prevax_sub_covidhospital_FALSE_ami <- cbind(
  method   = "lasso_X",
  name     = "cohort_prevax-sub_covidhospital_FALSE-ami",
  lasso_X_cohort_prevax_sub_covidhospital_FALSE_ami
)

lasso_X_cohort_prevax_sub_covidhospital_FALSE_ami <- subset(
  lasso_X_cohort_prevax_sub_covidhospital_FALSE_ami,
  select = -c(X)
)

lasso_X_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- cbind(
  method   = "lasso_X",
  name     = "cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs",
  lasso_X_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs
)

lasso_X_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- subset(
  lasso_X_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs,
  select = -c(X)
)

lasso_X_cohort_prevax_sub_covidhospital_TRUE_ami <- cbind(
  method   = "lasso_X",
  name     = "cohort_prevax-sub_covidhospital_TRUE-ami",
  lasso_X_cohort_prevax_sub_covidhospital_TRUE_ami
)

lasso_X_cohort_prevax_sub_covidhospital_TRUE_ami <- subset(
  lasso_X_cohort_prevax_sub_covidhospital_TRUE_ami,
  select = -c(X)
)

lasso_X_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs <- cbind(
  method   = "lasso_X",
  name     = "cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs",
  lasso_X_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)

lasso_X_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs <- subset(
  lasso_X_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs,
  select = -c(X)
)


# Stack lasso_X tables -----
print("Stack lasso_X tables")

lasso_X_all_cox_outputs <- rbind(
  lasso_X_cohort_prevax_main_ami,
  lasso_X_cohort_prevax_main_stroke_sahhs,
  lasso_X_cohort_prevax_sub_covidhospital_FALSE_ami,
  lasso_X_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs,
  lasso_X_cohort_prevax_sub_covidhospital_TRUE_ami,
  lasso_X_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Add method and name columns to each lasso_union cox table ----------------------
print("Add method and name columns to each lasso_union cox table")

lasso_union_cohort_prevax_main_ami <- cbind(
  method   = "lasso_union",
  name     = "cohort_prevax-main-ami",
  lasso_union_cohort_prevax_main_ami
)

lasso_union_cohort_prevax_main_ami <- subset(
  lasso_union_cohort_prevax_main_ami,
  select = -c(X)
)

lasso_union_cohort_prevax_main_stroke_sahhs <- cbind(
  method   = "lasso_union",
  name     = "cohort_prevax-main-stroke_sahhs",
  lasso_union_cohort_prevax_main_stroke_sahhs
)

lasso_union_cohort_prevax_main_stroke_sahhs <- subset(
  lasso_union_cohort_prevax_main_stroke_sahhs,
  select = -c(X)
)

lasso_union_cohort_prevax_sub_covidhospital_FALSE_ami <- cbind(
  method   = "lasso_union",
  name     = "cohort_prevax-sub_covidhospital_FALSE-ami",
  lasso_union_cohort_prevax_sub_covidhospital_FALSE_ami
)

lasso_union_cohort_prevax_sub_covidhospital_FALSE_ami <- subset(
  lasso_union_cohort_prevax_sub_covidhospital_FALSE_ami,
  select = -c(X)
)

lasso_union_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- cbind(
  method   = "lasso_union",
  name     = "cohort_prevax-sub_covidhospital_FALSE-stroke_sahhs",
  lasso_union_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs
)

lasso_union_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs <- subset(
  lasso_union_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs,
  select = -c(X)
)

lasso_union_cohort_prevax_sub_covidhospital_TRUE_ami <- cbind(
  method   = "lasso_union",
  name     = "cohort_prevax-sub_covidhospital_TRUE-ami",
  lasso_union_cohort_prevax_sub_covidhospital_TRUE_ami
)

lasso_union_cohort_prevax_sub_covidhospital_TRUE_ami <- subset(
  lasso_union_cohort_prevax_sub_covidhospital_TRUE_ami,
  select = -c(X)
)

lasso_union_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs <- cbind(
  method   = "lasso_union",
  name     = "cohort_prevax-sub_covidhospital_TRUE-stroke_sahhs",
  lasso_union_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)

lasso_union_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs <- subset(
  lasso_union_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs,
  select = -c(X)
)


# Stack lasso_union tables -----
print("Stack lasso_union tables")

lasso_union_all_cox_outputs <- rbind(
  lasso_union_cohort_prevax_main_ami,
  lasso_union_cohort_prevax_main_stroke_sahhs,
  lasso_union_cohort_prevax_sub_covidhospital_FALSE_ami,
  lasso_union_cohort_prevax_sub_covidhospital_FALSE_stroke_sahhs,
  lasso_union_cohort_prevax_sub_covidhospital_TRUE_ami,
  lasso_union_cohort_prevax_sub_covidhospital_TRUE_stroke_sahhs
)


# Save all tables -----
print("Save all tables")

all_cox_models_outputs <- rbind(
  fully_adjusted_all_cox_outputs,
  lasso_all_cox_outputs,
  lasso_X_all_cox_outputs,
  lasso_union_all_cox_outputs
)

write.csv(
  all_cox_models_outputs,
  paste0(makeout_dir, "all_cox_models_outputs.csv"),
  row.names = FALSE
)
