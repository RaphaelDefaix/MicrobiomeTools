#############################################################
## Run alpha diversity statistics
#############################################################

run_alpha_statistics <- function(
    alpha_long,
    group_var,
    day_var,
    alpha_indices
){

  message("Running alpha diversity statistics...")

  results <- list()

 days_to_compare <- levels(alpha_long[[day_var]])

days_to_compare <- days_to_compare[
  days_to_compare %in%
    unique(as.character(
      alpha_long[[day_var]][
        alpha_long[[group_var]] != "untreated"
      ]
    ))
]

  for(index in alpha_indices){

    message("-------------------------")
    message("Index : ", index)

    for(day in days_to_compare){

      df_test <- alpha_long %>%
        filter(
          Index == index,
          .data[[day_var]] == day,
          .data[[group_var]] != "untreated"
        ) %>%
        droplevels()

      message("   Day : ", day,
              " (n = ", nrow(df_test), ")")

    }

  }

  return(results)

}
