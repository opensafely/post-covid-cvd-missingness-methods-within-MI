preprocess_and_fit_model <- function(
  covariate_other,
  data
) {
  # Make binary variables logical ------------------------------------------------
  print("Make binary variables logical")

  var_bin <- colnames(data)[grepl("_bin_", colnames(data))]
  data[, var_bin] <- lapply(data[, var_bin], as.logical)

  # Make date variables dates ----------------------------------------------------
  print("Make date variables dates")

  var_date <- colnames(data)[grepl("_date", colnames(data))]
  data[, var_date] <- lapply(
    data[, var_date],
    function(x) as.Date(x, origin = "1970-01-01")
  )

  # Make categorical variables factors -------------------------------------------
  print("Make categorical variables factors")

  var_cat <- colnames(data)[grepl("_cat_", colnames(data))]
  data[, var_cat] <- lapply(data[, var_cat], as.factor)

  # Make numerical variables numerical -------------------------------------------
  print(" Make numerical variables numerical")

  var_num <- colnames(data)[grepl("_num_", colnames(data))]
  data[, var_num] <- lapply(data[, var_num], as.numeric)

  # Study date variables ---------------------------------------------------------
  print("Study date variables")

  data <- dplyr::rename(
    data,
    "outcome" = tidyselect::all_of("out_date"),
    "exposure" = tidyselect::all_of("exp_date")
  )

  cox_start <- "index_date"
  cox_stop  <- "end_date_outcome"

  study_start <- "2020-01-01"
  study_stop  <- "2024-04-30"

  # cox_start <- gsub(out_date, "outcome", cox_start)
  # cox_start <- gsub(exp_date, "exposure", cox_start)

  # cox_stop <- gsub(out_date, "outcome", cox_stop)
  # cox_stop <- gsub(exp_date, "exposure", cox_stop)

  # data$cox_start <- as.Date(cox_start)
  # data$cox_stop  <- as.Date(cox_stop)

  data$study_start <- as.Date(study_start)
  data$study_stop <- as.Date(study_stop)

  data$fup_start <- do.call(
    pmax,
    c(data[, c("study_start", cox_start)], list(na.rm = TRUE))
  )

  data$fup_stop <- do.call(
    pmin,
    c(data[, c("study_stop", cox_stop, "outcome")], list(na.rm = TRUE))
  )

  data <- data[data$fup_stop >= data$fup_start, ]

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

  # Remove exposures and outcomes outside follow-up ------------------------------
  print("Remove exposures and outcomes outside follow-up")

  print(paste0(
    "Exposure data range: ",
    min(data$exposure, na.rm = TRUE),
    " to ",
    max(data$exposure, na.rm = TRUE)
  ))
  print(paste0(
    "Outcome data range: ",
    min(data$outcome, na.rm = TRUE),
    " to ",
    max(data$outcome, na.rm = TRUE)
  ))

  data <- data %>%
    dplyr::mutate(
      exposure = replace(
        exposure,
        which(exposure > fup_stop | exposure < fup_start),
        NA
      ),
      outcome = replace(
        outcome,
        which(outcome > fup_stop | outcome < fup_start),
        NA
      )
    )

  print(paste0(
    "Exposure data range: ",
    min(data$exposure, na.rm = TRUE),
    " to ",
    max(data$exposure, na.rm = TRUE)
  ))
  print(paste0(
    "Outcome data range: ",
    min(data$outcome, na.rm = TRUE),
    " to ",
    max(data$outcome, na.rm = TRUE)
  ))

  # Make indicator variable for outcome status -----------------------------------
  print("Make indicator variable for outcome status")

  data$outcome_status <- data$outcome == data$fup_stop &
    !is.na(data$outcome) &
    !is.na(data$fup_stop)

  print(table(data$outcome_status))

  data_surv <- survival_data_setup(
    data             = data,
    cut_points     = cut_points,
    episode_labels = episode_labels
  )

  episode_info <- get_episode_info(
    data             = data_surv,
    cut_points     = cut_points,
    episode_labels = episode_labels,
    ipw            = FALSE
  )

  # MODEL CODE STARTS HERE
  data                  <- data_surv
  time_periods        <- episode_info[
    episode_info$time_period != "days_pre",
  ]$time_period
  covariates          <- covariate_other
  strata              <- "strat_cat_region"
  age_spline          <- TRUE
  covariate_removed   <- NULL
  covariate_collapsed <- NULL
  ipw                 <- FALSE

  # Define model formula -------------------------------------------------------
  print("Define model formula")

  surv_formula <- paste0(
    "Surv(tstart, tstop, outcome_status) ~ ",
    paste(time_periods, collapse = " + "),
    ifelse("cov_cat_sex" %in% colnames(data), " + cov_cat_sex", ""),
    ifelse(
      is.null(strata),
      "",
      paste(" +", paste0("rms::strat(", strata, ")"), collapse = " + ")
    ),
    ifelse(isTRUE(ipw), " + cluster(patient_id)", "")
  )

  # Add age covariate, specifying knot placement for age spline if applicable --

  if ("cov_num_age" %in% colnames(data)) {
    print("Add age covariate")

    if (age_spline == TRUE) {
      print("Specify knot placement for age spline")

      knot_placement <- as.numeric(quantile(
        data$cov_num_age,
        probs = c(0.1, 0.5, 0.9)
      ))

      print(paste0(
        "Knots will be placed at: ",
        paste0(knot_placement, collapse = ", ")
      ))

      surv_formula <- paste0(
        surv_formula,
        " + rms::rcs(cov_num_age, parms=knot_placement)"
      )
    } else {
      surv_formula <- paste0(surv_formula, " + cov_num_age + cov_num_age_sq")
    }
  }

  print(surv_formula)

  # Fit Cox model ----------------------------------------------------------------
  print("Fit Cox model")

  dd <<- rms::datadist(data)

  withr::local_options(list(
    datadist = "dd",
    contrasts = c("contr.treatment", "contr.treatment")
  ))

  if (ipw == TRUE) {
    N_obs_in <- nrow(data)

    fit_cox_model <- rms::cph(
      formula = as.formula(surv_formula),
      data = data,
      weight = data$cox_weight,
      method = "breslow",
      surv = TRUE,
      x = TRUE,
      y = TRUE
    )

    N_obs_out <- sum(fit_cox_model$n)
  } else {
    N_obs_in <- nrow(data)

    fit_cox_model <- rms::cph(
      formula = as.formula(surv_formula),
      data = data,
      method = "breslow",
      surv = TRUE,
      x = TRUE,
      y = TRUE
    )

    N_obs_out <- sum(fit_cox_model$n)
  }

  print(fit_cox_model)

  # Format results ---------------------------------------------------------------
  print("Format results")

  results <- data.frame(
    term = names(fit_cox_model$coefficients),
    lnhr = fit_cox_model$coefficients,
    se_lnhr = sqrt(diag(vcov(fit_cox_model))),
    model = "mdl_age_sex",
    surv_formula = surv_formula,
    covariate_removed = "",
    covariate_collapsed = "",
    obs_warning = ifelse(
      N_obs_in == N_obs_out,
      "",
      paste0(
        N_obs_in,
        " observations provided. ",
        N_obs_out,
        " observations used."
      )
    ),
    stringsAsFactors = FALSE
  )

  row.names(results) <- NULL

  # If covariates are specified, run an additional model including them --------

  covariates <- setdiff(covariates, covariate_removed)

  if (!is.null(covariates) && length(covariates) > 0) {
    # Add covariates to model formula ------------------------------------------
    print("Add covariates to model formula")

    surv_formula_adj <- paste0(
      surv_formula,
      " + ",
      paste(covariates, collapse = " + ")
    )

    print(surv_formula_adj)

    # Fit Cox model ------------------------------------------------------------
    print("Fit Cox model with covariates")

    dd_adj <<- rms::datadist(data)

    withr::local_options(list(
      datadist = "dd_adj",
      contrasts = c("contr.treatment", "contr.treatment")
    ))

    if (ipw == TRUE) {
      N_obs_in <- nrow(data)

      fit_cox_model_adj <- rms::cph(
        formula = as.formula(surv_formula_adj),
        data = data,
        weight = data$cox_weight,
        method = "breslow",
        surv = TRUE,
        x = TRUE,
        y = TRUE
      )

      N_obs_out <- sum(fit_cox_model$n)
    } else {
      N_obs_in <- nrow(data)

      fit_cox_model_adj <- rms::cph(
        formula = as.formula(surv_formula_adj),
        data = data,
        method = "breslow",
        surv = TRUE,
        x = TRUE,
        y = TRUE
      )

      N_obs_out <- sum(fit_cox_model$n)
    }

    print(fit_cox_model_adj)

    # Format results ---------------------------------------------------------------
    print("Format results")

    results_adj <- data.frame(
      term = names(fit_cox_model_adj$coefficients),
      lnhr = fit_cox_model_adj$coefficients,
      se_lnhr = sqrt(diag(vcov(fit_cox_model_adj))),
      model = "mdl_max_adj",
      surv_formula = surv_formula_adj,
      covariate_removed = paste0(covariate_removed, collapse = ";"),
      covariate_collapsed = paste0(covariate_collapsed, collapse = ";"),
      obs_warning = ifelse(
        N_obs_in == N_obs_out,
        "",
        paste0(
          N_obs_in,
          " observations provided. ",
          N_obs_out,
          " observations used."
        )
      ),
      stringsAsFactors = FALSE
    )

    row.names(results_adj) <- NULL

    # Bind to other results ----------------------------------------------------
    print("Bind to other results")

    results <- rbind(results, results_adj)
  }

  # Return results -------------------------------------------------------------
  print("Return results")

  return(results)
  print(summary(results))
}
