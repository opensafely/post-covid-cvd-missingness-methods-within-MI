get_cox_results <- function(list_of_cox_models = NULL, surv_formula = NULL) {
  list_of_cox_results <- list()

  for (i in 1:get_number_of_imputed_datasets()) {
    model_i   <- list_of_cox_models[[i]]

    results_i <- data.frame(
      term                = names(model_i$coefficients),
      lnhr                = model_i$coefficients,
      se_lnhr             = sqrt(diag(vcov(model_i))),
      model               = "mdl_max_adj",
      surv_formula        = surv_formula,
      covariate_removed   = "",
      covariate_collapsed = "",
      stringsAsFactors    = FALSE
    )
    rownames(results_i) <- NULL # duplicate term names

    list_of_cox_results[[i]] <- results_i
  }

  return(list_of_cox_results)
}
