print_out_cox_summaries <- function(list_of_cox_results = NULL) {
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

  print(dataset_1)
  print(dataset_2)
  print(dataset_3)
  print(dataset_4)
  print(dataset_5)
  print(dataset_6)
  print(dataset_7)
  print(dataset_8)
  print(dataset_9)
  print(dataset_10)
  stop("m=10 datasets")

  return (NULL)
}
