manual_pool_RR <- function(list_of_cox_results = NULL) {
  library(dplyr)

  day_terms <- c(
    "days0_1", "days1_28", "days28_196", "days196_364", "days364_714", "days714_1582"
  )

  # Extract results per dataset ---------------------------------------------
  print("Extract results per dataset")

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
  

  # Mean lnhr from each dataset -------------------------------
  print("Mean lnhr from each dataset")

  pooled_lnhr_days0_1 <- mean(
    (dataset_1[which(dataset_1$term == "days0_1"), "lnhr"]),
    (dataset_2[which(dataset_2$term == "days0_1"), "lnhr"]),
    (dataset_3[which(dataset_3$term == "days0_1"), "lnhr"]),
    (dataset_4[which(dataset_4$term == "days0_1"), "lnhr"]),
    (dataset_5[which(dataset_5$term == "days0_1"), "lnhr"]),
    (dataset_6[which(dataset_6$term == "days0_1"), "lnhr"]),
    (dataset_7[which(dataset_7$term == "days0_1"), "lnhr"]),
    (dataset_8[which(dataset_8$term == "days0_1"), "lnhr"]),
    (dataset_9[which(dataset_9$term == "days0_1"), "lnhr"]),
    (dataset_10[which(dataset_10$term == "days0_1"), "lnhr"])
  )

  pooled_lnhr_days1_28 <- mean(
    (dataset_1[which(dataset_1$term == "days1_28"), "lnhr"]),
    (dataset_2[which(dataset_2$term == "days1_28"), "lnhr"]),
    (dataset_3[which(dataset_3$term == "days1_28"), "lnhr"]),
    (dataset_4[which(dataset_4$term == "days1_28"), "lnhr"]),
    (dataset_5[which(dataset_5$term == "days1_28"), "lnhr"]),
    (dataset_6[which(dataset_6$term == "days1_28"), "lnhr"]),
    (dataset_7[which(dataset_7$term == "days1_28"), "lnhr"]),
    (dataset_8[which(dataset_8$term == "days1_28"), "lnhr"]),
    (dataset_9[which(dataset_9$term == "days1_28"), "lnhr"]),
    (dataset_10[which(dataset_10$term == "days1_28"), "lnhr"])
  )

  pooled_lnhr_days28_196 <- mean(
    (dataset_1[which(dataset_1$term == "days28_196"), "lnhr"]),
    (dataset_2[which(dataset_2$term == "days28_196"), "lnhr"]),
    (dataset_3[which(dataset_3$term == "days28_196"), "lnhr"]),
    (dataset_4[which(dataset_4$term == "days28_196"), "lnhr"]),
    (dataset_5[which(dataset_5$term == "days28_196"), "lnhr"]),
    (dataset_6[which(dataset_6$term == "days28_196"), "lnhr"]),
    (dataset_7[which(dataset_7$term == "days28_196"), "lnhr"]),
    (dataset_8[which(dataset_8$term == "days28_196"), "lnhr"]),
    (dataset_9[which(dataset_9$term == "days28_196"), "lnhr"]),
    (dataset_10[which(dataset_10$term == "days28_196"), "lnhr"])
  )

  pooled_lnhr_days196_364 <- mean(
    (dataset_1[which(dataset_1$term == "days196_364"), "lnhr"]),
    (dataset_2[which(dataset_2$term == "days196_364"), "lnhr"]),
    (dataset_3[which(dataset_3$term == "days196_364"), "lnhr"]),
    (dataset_4[which(dataset_4$term == "days196_364"), "lnhr"]),
    (dataset_5[which(dataset_5$term == "days196_364"), "lnhr"]),
    (dataset_6[which(dataset_6$term == "days196_364"), "lnhr"]),
    (dataset_7[which(dataset_7$term == "days196_364"), "lnhr"]),
    (dataset_8[which(dataset_8$term == "days196_364"), "lnhr"]),
    (dataset_9[which(dataset_9$term == "days196_364"), "lnhr"]),
    (dataset_10[which(dataset_10$term == "days196_364"), "lnhr"])
  )

  pooled_lnhr_days364_714 <- mean(
    (dataset_1[which(dataset_1$term == "days364_714"), "lnhr"]),
    (dataset_2[which(dataset_2$term == "days364_714"), "lnhr"]),
    (dataset_3[which(dataset_3$term == "days364_714"), "lnhr"]),
    (dataset_4[which(dataset_4$term == "days364_714"), "lnhr"]),
    (dataset_5[which(dataset_5$term == "days364_714"), "lnhr"]),
    (dataset_6[which(dataset_6$term == "days364_714"), "lnhr"]),
    (dataset_7[which(dataset_7$term == "days364_714"), "lnhr"]),
    (dataset_8[which(dataset_8$term == "days364_714"), "lnhr"]),
    (dataset_9[which(dataset_9$term == "days364_714"), "lnhr"]),
    (dataset_10[which(dataset_10$term == "days364_714"), "lnhr"])
  )

  pooled_lnhr_days714_1582 <- mean(
    (dataset_1[which(dataset_1$term == "days714_1582"), "lnhr"]),
    (dataset_2[which(dataset_2$term == "days714_1582"), "lnhr"]),
    (dataset_3[which(dataset_3$term == "days714_1582"), "lnhr"]),
    (dataset_4[which(dataset_4$term == "days714_1582"), "lnhr"]),
    (dataset_5[which(dataset_5$term == "days714_1582"), "lnhr"]),
    (dataset_6[which(dataset_6$term == "days714_1582"), "lnhr"]),
    (dataset_7[which(dataset_7$term == "days714_1582"), "lnhr"]),
    (dataset_8[which(dataset_8$term == "days714_1582"), "lnhr"]),
    (dataset_9[which(dataset_9$term == "days714_1582"), "lnhr"]),
    (dataset_10[which(dataset_10$term == "days714_1582"), "lnhr"])
  )


  # Within variance from each dataset -------------------------------
  print("Within variance from each dataset")

  within_variance_days0_1      <- mean(
    (dataset_1[which(dataset_1$term == "days0_1"), "se_lnhr"])^2,
    (dataset_2[which(dataset_2$term == "days0_1"), "se_lnhr"])^2,
    (dataset_3[which(dataset_3$term == "days0_1"), "se_lnhr"])^2,
    (dataset_4[which(dataset_4$term == "days0_1"), "se_lnhr"])^2,
    (dataset_5[which(dataset_5$term == "days0_1"), "se_lnhr"])^2,
    (dataset_6[which(dataset_6$term == "days0_1"), "se_lnhr"])^2,
    (dataset_7[which(dataset_7$term == "days0_1"), "se_lnhr"])^2,
    (dataset_8[which(dataset_8$term == "days0_1"), "se_lnhr"])^2,
    (dataset_9[which(dataset_9$term == "days0_1"), "se_lnhr"])^2,
    (dataset_10[which(dataset_10$term == "days0_1"), "se_lnhr"])^2
  )
  within_variance_days1_28      <- mean(
    (dataset_1[which(dataset_1$term == "days1_28"), "se_lnhr"])^2,
    (dataset_2[which(dataset_2$term == "days1_28"), "se_lnhr"])^2,
    (dataset_3[which(dataset_3$term == "days1_28"), "se_lnhr"])^2,
    (dataset_4[which(dataset_4$term == "days1_28"), "se_lnhr"])^2,
    (dataset_5[which(dataset_5$term == "days1_28"), "se_lnhr"])^2,
    (dataset_6[which(dataset_6$term == "days1_28"), "se_lnhr"])^2,
    (dataset_7[which(dataset_7$term == "days1_28"), "se_lnhr"])^2,
    (dataset_8[which(dataset_8$term == "days1_28"), "se_lnhr"])^2,
    (dataset_9[which(dataset_9$term == "days1_28"), "se_lnhr"])^2,
    (dataset_10[which(dataset_10$term == "days1_28"), "se_lnhr"])^2
  )
  within_variance_days28_196      <- mean(
    (dataset_1[which(dataset_1$term == "days28_196"), "se_lnhr"])^2,
    (dataset_2[which(dataset_2$term == "days28_196"), "se_lnhr"])^2,
    (dataset_3[which(dataset_3$term == "days28_196"), "se_lnhr"])^2,
    (dataset_4[which(dataset_4$term == "days28_196"), "se_lnhr"])^2,
    (dataset_5[which(dataset_5$term == "days28_196"), "se_lnhr"])^2,
    (dataset_6[which(dataset_6$term == "days28_196"), "se_lnhr"])^2,
    (dataset_7[which(dataset_7$term == "days28_196"), "se_lnhr"])^2,
    (dataset_8[which(dataset_8$term == "days28_196"), "se_lnhr"])^2,
    (dataset_9[which(dataset_9$term == "days28_196"), "se_lnhr"])^2,
    (dataset_10[which(dataset_10$term == "days28_196"), "se_lnhr"])^2
  )
  within_variance_days196_364      <- mean(
    (dataset_1[which(dataset_1$term == "days196_364"), "se_lnhr"])^2,
    (dataset_2[which(dataset_2$term == "days196_364"), "se_lnhr"])^2,
    (dataset_3[which(dataset_3$term == "days196_364"), "se_lnhr"])^2,
    (dataset_4[which(dataset_4$term == "days196_364"), "se_lnhr"])^2,
    (dataset_5[which(dataset_5$term == "days196_364"), "se_lnhr"])^2,
    (dataset_6[which(dataset_6$term == "days196_364"), "se_lnhr"])^2,
    (dataset_7[which(dataset_7$term == "days196_364"), "se_lnhr"])^2,
    (dataset_8[which(dataset_8$term == "days196_364"), "se_lnhr"])^2,
    (dataset_9[which(dataset_9$term == "days196_364"), "se_lnhr"])^2,
    (dataset_10[which(dataset_10$term == "days196_364"), "se_lnhr"])^2
  )
  within_variance_days364_714      <- mean(
    (dataset_1[which(dataset_1$term == "days364_714"), "se_lnhr"])^2,
    (dataset_2[which(dataset_2$term == "days364_714"), "se_lnhr"])^2,
    (dataset_3[which(dataset_3$term == "days364_714"), "se_lnhr"])^2,
    (dataset_4[which(dataset_4$term == "days364_714"), "se_lnhr"])^2,
    (dataset_5[which(dataset_5$term == "days364_714"), "se_lnhr"])^2,
    (dataset_6[which(dataset_6$term == "days364_714"), "se_lnhr"])^2,
    (dataset_7[which(dataset_7$term == "days364_714"), "se_lnhr"])^2,
    (dataset_8[which(dataset_8$term == "days364_714"), "se_lnhr"])^2,
    (dataset_9[which(dataset_9$term == "days364_714"), "se_lnhr"])^2,
    (dataset_10[which(dataset_10$term == "days364_714"), "se_lnhr"])^2
  )
  within_variance_days714_1582      <- mean(
    (dataset_1[which(dataset_1$term == "days714_1582"), "se_lnhr"])^2,
    (dataset_2[which(dataset_2$term == "days714_1582"), "se_lnhr"])^2,
    (dataset_3[which(dataset_3$term == "days714_1582"), "se_lnhr"])^2,
    (dataset_4[which(dataset_4$term == "days714_1582"), "se_lnhr"])^2,
    (dataset_5[which(dataset_5$term == "days714_1582"), "se_lnhr"])^2,
    (dataset_6[which(dataset_6$term == "days714_1582"), "se_lnhr"])^2,
    (dataset_7[which(dataset_7$term == "days714_1582"), "se_lnhr"])^2,
    (dataset_8[which(dataset_8$term == "days714_1582"), "se_lnhr"])^2,
    (dataset_9[which(dataset_9$term == "days714_1582"), "se_lnhr"])^2,
    (dataset_10[which(dataset_10$term == "days714_1582"), "se_lnhr"])^2
  )


  # Between variance from each dataset -------------------------------
  print("Between variance from each dataset")

  between_variance_days0_1      <- var(c(
    (dataset_1[which(dataset_1$term == "days0_1"), "lnhr"]),
    (dataset_2[which(dataset_2$term == "days0_1"), "lnhr"]),
    (dataset_3[which(dataset_3$term == "days0_1"), "lnhr"]),
    (dataset_4[which(dataset_4$term == "days0_1"), "lnhr"]),
    (dataset_5[which(dataset_5$term == "days0_1"), "lnhr"]),
    (dataset_6[which(dataset_6$term == "days0_1"), "lnhr"]),
    (dataset_7[which(dataset_7$term == "days0_1"), "lnhr"]),
    (dataset_8[which(dataset_8$term == "days0_1"), "lnhr"]),
    (dataset_9[which(dataset_9$term == "days0_1"), "lnhr"]),
    (dataset_10[which(dataset_10$term == "days0_1"), "lnhr"])
  ))
  between_variance_days1_28      <- var(c(
    (dataset_1[which(dataset_1$term == "days1_28"), "lnhr"]),
    (dataset_2[which(dataset_2$term == "days1_28"), "lnhr"]),
    (dataset_3[which(dataset_3$term == "days1_28"), "lnhr"]),
    (dataset_4[which(dataset_4$term == "days1_28"), "lnhr"]),
    (dataset_5[which(dataset_5$term == "days1_28"), "lnhr"]),
    (dataset_6[which(dataset_6$term == "days1_28"), "lnhr"]),
    (dataset_7[which(dataset_7$term == "days1_28"), "lnhr"]),
    (dataset_8[which(dataset_8$term == "days1_28"), "lnhr"]),
    (dataset_9[which(dataset_9$term == "days1_28"), "lnhr"]),
    (dataset_10[which(dataset_10$term == "days1_28"), "lnhr"])
  ))
  between_variance_days28_196      <- var(c(
    (dataset_1[which(dataset_1$term == "days28_196"), "lnhr"]),
    (dataset_2[which(dataset_2$term == "days28_196"), "lnhr"]),
    (dataset_3[which(dataset_3$term == "days28_196"), "lnhr"]),
    (dataset_4[which(dataset_4$term == "days28_196"), "lnhr"]),
    (dataset_5[which(dataset_5$term == "days28_196"), "lnhr"]),
    (dataset_6[which(dataset_6$term == "days28_196"), "lnhr"]),
    (dataset_7[which(dataset_7$term == "days28_196"), "lnhr"]),
    (dataset_8[which(dataset_8$term == "days28_196"), "lnhr"]),
    (dataset_9[which(dataset_9$term == "days28_196"), "lnhr"]),
    (dataset_10[which(dataset_10$term == "days28_196"), "lnhr"])
  ))
  between_variance_days196_364      <- var(c(
    (dataset_1[which(dataset_1$term == "days196_364"), "lnhr"]),
    (dataset_2[which(dataset_2$term == "days196_364"), "lnhr"]),
    (dataset_3[which(dataset_3$term == "days196_364"), "lnhr"]),
    (dataset_4[which(dataset_4$term == "days196_364"), "lnhr"]),
    (dataset_5[which(dataset_5$term == "days196_364"), "lnhr"]),
    (dataset_6[which(dataset_6$term == "days196_364"), "lnhr"]),
    (dataset_7[which(dataset_7$term == "days196_364"), "lnhr"]),
    (dataset_8[which(dataset_8$term == "days196_364"), "lnhr"]),
    (dataset_9[which(dataset_9$term == "days196_364"), "lnhr"]),
    (dataset_10[which(dataset_10$term == "days196_364"), "lnhr"])
  ))
  between_variance_days364_714      <- var(c(
    (dataset_1[which(dataset_1$term == "days364_714"), "lnhr"]),
    (dataset_2[which(dataset_2$term == "days364_714"), "lnhr"]),
    (dataset_3[which(dataset_3$term == "days364_714"), "lnhr"]),
    (dataset_4[which(dataset_4$term == "days364_714"), "lnhr"]),
    (dataset_5[which(dataset_5$term == "days364_714"), "lnhr"]),
    (dataset_6[which(dataset_6$term == "days364_714"), "lnhr"]),
    (dataset_7[which(dataset_7$term == "days364_714"), "lnhr"]),
    (dataset_8[which(dataset_8$term == "days364_714"), "lnhr"]),
    (dataset_9[which(dataset_9$term == "days364_714"), "lnhr"]),
    (dataset_10[which(dataset_10$term == "days364_714"), "lnhr"])
  ))
  between_variance_days714_1582      <- var(c(
    (dataset_1[which(dataset_1$term == "days714_1582"), "lnhr"]),
    (dataset_2[which(dataset_2$term == "days714_1582"), "lnhr"]),
    (dataset_3[which(dataset_3$term == "days714_1582"), "lnhr"]),
    (dataset_4[which(dataset_4$term == "days714_1582"), "lnhr"]),
    (dataset_5[which(dataset_5$term == "days714_1582"), "lnhr"]),
    (dataset_6[which(dataset_6$term == "days714_1582"), "lnhr"]),
    (dataset_7[which(dataset_7$term == "days714_1582"), "lnhr"]),
    (dataset_8[which(dataset_8$term == "days714_1582"), "lnhr"]),
    (dataset_9[which(dataset_9$term == "days714_1582"), "lnhr"]),
    (dataset_10[which(dataset_10$term == "days714_1582"), "lnhr"])
  ))


  # Total variance from across all datasets -------------------------------
  print("Total variance from across all datasets")

  all_variance_days0_1      <- within_variance_days0_1 + between_variance_days0_1 + (between_variance_days0_1 / 10)
  all_variance_days1_28     <- within_variance_days1_28 + between_variance_days1_28 + (between_variance_days1_28 / 10)
  all_variance_days28_196   <- within_variance_days28_196 + between_variance_days28_196 + (between_variance_days28_196 / 10)
  all_variance_days196_364  <- within_variance_days196_364 + between_variance_days196_364 + (between_variance_days196_364 / 10)
  all_variance_days364_714  <- within_variance_days364_714 + between_variance_days364_714 + (between_variance_days364_714 / 10)
  all_variance_days714_1582 <- within_variance_days714_1582 + between_variance_days714_1582 + (between_variance_days714_1582 / 10)


  # Total SE from across all datasets -------------------------------
  print("Total SE from across all datasets")

  all_SE_days0_1      <- sqrt(all_variance_days0_1)
  all_SE_days1_28     <- sqrt(all_variance_days1_28)
  all_SE_days28_196   <- sqrt(all_variance_days28_196)
  all_SE_days196_364  <- sqrt(all_variance_days196_364)
  all_SE_days364_714  <- sqrt(all_variance_days364_714)
  all_SE_days714_1582 <- sqrt(all_variance_days714_1582)

  pooled_model_results <- data.frame(
    term = c(
      "days0_1",
      "days1_28",
      "days28_196",
      "days196_364",
      "days364_714",
      "days714_1582"
    ),
    lnhr = c(
      pooled_lnhr_days0_1,
      pooled_lnhr_days1_28,
      pooled_lnhr_days28_196,
      pooled_lnhr_days196_364,
      pooled_lnhr_days364_714,
      pooled_lnhr_days714_1582
    ),
    se_lnhr = c(
      all_SE_days0_1,
      all_SE_days1_28,
      all_SE_days28_196,
      all_SE_days196_364,
      all_SE_days364_714,
      all_SE_days714_1582
    )
  )

  return(pooled_model_results)
}
