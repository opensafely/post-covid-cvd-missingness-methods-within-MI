manual_pool_RR <- function(list_of_cox_results = NULL) {
  # extract data from first set of results
  columns <- colnames(list_of_cox_results[[1]])
  term    <- (list_of_cox_results[[1]])$term
  model   <- (list_of_cox_results[[1]])$model
  surv_formula          <- (list_of_cox_results[[1]])$surv_formula
  covariate_removed     <- (list_of_cox_results[[1]])$covariate_removed
  covariate_collapsed   <- (list_of_cox_results[[1]])$covariate_collapsed

  # initialise dataframe of pooled results
  pooled_model_results <- data.frame(
    matrix(data = NA, nrow = length(term), ncol = length(columns))
  )
  colnames(pooled_model_results) <- columns

  # set any terms that are constant
  pooled_model_results$term  <- term
  pooled_model_results$model <- model
  pooled_model_results$surv_formula        <- surv_formula
  pooled_model_results$covariate_removed   <- covariate_removed
  pooled_model_results$covariate_collapsed <- covariate_collapsed

  # initialise dataframe of coefs
  all_model_lnhr <- data.frame(
    matrix(nrow = length(term), ncol = get_number_of_imputed_datasets())
  )
  rownames(all_model_lnhr) <- term

  # initialise dataframe of SEs
  all_model_se_lnhr <- data.frame(
    matrix(nrow = length(term), ncol = get_number_of_imputed_datasets())
  )
  rownames(all_model_se_lnhr) <- term

  # extract data from listy of dataframes
  for (i in c(1:get_number_of_imputed_datasets())) {
    all_model_lnhr[, i]    <- (list_of_cox_results[[i]])$lnhr
    all_model_se_lnhr[, i] <- (list_of_cox_results[[i]])$se_lnhr
  }

  # within imputation variance
  all_model_within_variance   <- all_model_se_lnhr^2 # element-wise operation
  mean_estimate               <- rowMeans(all_model_lnhr)
  
  # deviance of each lnhr estimate from mean estimate
  deviance <- data.frame(
    matrix(nrow = length(term), ncol = get_number_of_imputed_datasets())
  )
  rownames(deviance) <- term
  
  for (i in c(1:length(term))) {
    for (j in c(1:get_number_of_imputed_datasets())) {
      deviance[i, j] <- (all_model_lnhr[i, j] - mean_estimate[j])
    }
  }
  
  # between and total variances
  sq_deviance                 <- (deviance)^2
  all_model_between_variance  <- ((rowSums(sq_deviance)) / (get_number_of_imputed_datasets() - 1))
  all_model_total_variance    <- all_model_within_variance + (all_model_between_variance) + (all_model_between_variance / get_number_of_imputed_datasets())
  
  # extract final values
  pooled_lnhr     <- rowMeans(all_model_lnhr)
  pooled_variance <- rowMeans(all_model_total_variance)
  pooled_lnhr_se  <- sqrt(pooled_variance) # element-wise operation

  # set variable columns (coefs and SEs)
  pooled_model_results$lnhr    <- unname(pooled_lnhr)
  pooled_model_results$se_lnhr <- unname(pooled_lnhr_se)

  return(pooled_model_results)
}
