manual_pool_RR <- function(list_of_cox_results = NULL) {
  library(dplyr)

  day_terms <- c(
    "days0_1", "days1_28", "days28_196", "days196_364", "days364_714", "days714_1582"
  )

  dataset_1 <- list_of_cox_results[[1]] %>%
    dplyr::select(term, lnhr, se_lnhr) %>%
    dplyr::filter(term %in% day_terms)
  
  dataset_2 <- list_of_cox_results[[2]] %>%
    dplyr::select(term, lnhr, se_lnhr) %>%
    dplyr::filter(term %in% day_terms)

  dataset_3 <- list_of_cox_results[[3]] %>%
    dplyr::select(term, lnhr, se_lnhr) %>%
    dplyr::filter(term %in% day_terms)
  
  dataset_4 <- list_of_cox_results[[4]] %>%
    dplyr::select(term, lnhr, se_lnhr) %>%
    dplyr::filter(term %in% day_terms)
  
  dataset_5 <- list_of_cox_results[[5]] %>%
    dplyr::select(term, lnhr, se_lnhr) %>%
    dplyr::filter(term %in% day_terms)
  
  dataset_6 <- list_of_cox_results[[6]] %>%
    dplyr::select(term, lnhr, se_lnhr) %>%
    dplyr::filter(term %in% day_terms)
  
  dataset_7 <- list_of_cox_results[[7]] %>%
    dplyr::select(term, lnhr, se_lnhr) %>%
    dplyr::filter(term %in% day_terms)
  
  dataset_8 <- list_of_cox_results[[8]] %>%
    dplyr::select(term, lnhr, se_lnhr) %>%
    dplyr::filter(term %in% day_terms)
  
  dataset_9 <- list_of_cox_results[[9]] %>%
    dplyr::select(term, lnhr, se_lnhr) %>%
    dplyr::filter(term %in% day_terms)
  
  dataset_10 <- list_of_cox_results[[10]] %>%
    dplyr::select(term, lnhr, se_lnhr) %>%
    dplyr::filter(term %in% day_terms)
  
  within_variance_days0_1      <- mean(
    (dataset_1["days1_28", se_lnhr])^2 +
    (dataset_2["days1_28", se_lnhr])^2
  )
  print(within_variance_days0_1)
  stop("!!!")

  within_variance_days1_28     <- NULL
  within_variance_days28_196   <- NULL
  within_variance_days196_364  <- NULL
  within_variance_days364_714  <- NULL
  within_variance_days714_1582 <- NULL

  between_variance_days0_1      <- NULL
  between_variance_days1_28     <- NULL
  between_variance_days28_196   <- NULL
  between_variance_days196_364  <- NULL
  between_variance_days364_714  <- NULL
  between_variance_days714_1582 <- NULL

  all_variance_days0_1      <- NULL
  all_variance_days1_28     <- NULL
  all_variance_days28_196   <- NULL
  all_variance_days196_364  <- NULL
  all_variance_days364_714  <- NULL
  all_variance_days714_1582 <- NULL

  return(pooled_model_results)
}
