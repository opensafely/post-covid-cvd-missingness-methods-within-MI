preprocess_surv_data <- function(
  df
) {
  # Make binary variables logical ------------------------------------------------
  print("Make binary variables logical")

  var_bin <- colnames(df)[grepl("_bin_", colnames(df))]
  df[, var_bin] <- lapply(df[, var_bin], as.logical)

  # Make date variables dates ----------------------------------------------------
  print("Make date variables dates")

  var_date <- colnames(df)[grepl("_date", colnames(df))]
  df[, var_date] <- lapply(
    df[, var_date],
    function(x) as.Date(x, origin = "1970-01-01")
  )

  # Make categorical variables factors -------------------------------------------
  print("Make categorical variables factors")

  var_cat <- colnames(df)[grepl("_cat_", colnames(df))]
  df[, var_cat] <- lapply(df[, var_cat], as.factor)

  # Make numerical variables numerical -------------------------------------------
  print(" Make numerical variables numerical")

  var_num <- colnames(df)[grepl("_num_", colnames(df))]
  df[, var_num] <- lapply(df[, var_num], as.numeric)

  # Study date variables ---------------------------------------------------------
  print("Study date variables")

  df <- dplyr::rename(
    df,
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

  # df$cox_start <- as.Date(cox_start)
  # df$cox_stop  <- as.Date(cox_stop)

  df$study_start <- as.Date(study_start)
  df$study_stop <- as.Date(study_stop)

  df$fup_start <- do.call(
    pmax,
    c(df[, c("study_start", cox_start)], list(na.rm = TRUE))
  )

  df$fup_stop <- do.call(
    pmin,
    c(df[, c("study_stop", cox_stop, "outcome")], list(na.rm = TRUE))
  )

  df <- df[df$fup_stop >= df$fup_start, ]

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
    min(df$exposure, na.rm = TRUE),
    " to ",
    max(df$exposure, na.rm = TRUE)
  ))
  print(paste0(
    "Outcome data range: ",
    min(df$outcome, na.rm = TRUE),
    " to ",
    max(df$outcome, na.rm = TRUE)
  ))

  df <- df %>%
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
    min(df$exposure, na.rm = TRUE),
    " to ",
    max(df$exposure, na.rm = TRUE)
  ))
  print(paste0(
    "Outcome data range: ",
    min(df$outcome, na.rm = TRUE),
    " to ",
    max(df$outcome, na.rm = TRUE)
  ))

  # Make indicator variable for outcome status -----------------------------------
  print("Make indicator variable for outcome status")

  df$outcome_status <- df$outcome == df$fup_stop &
    !is.na(df$outcome) &
    !is.na(df$fup_stop)

  print(table(df$outcome_status))

  data_surv <- survival_data_setup(
    df             = df,
    cut_points     = cut_points,
    episode_labels = episode_labels
  )

  episode_info <- get_episode_info(
    df             = data_surv,
    cut_points     = cut_points,
    episode_labels = episode_labels,
    ipw            = FALSE
  )

  return(list(
    data_surv    = data_surv,
    episode_info = episode_info
  ))
}
