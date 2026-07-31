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

  ## Days with at least two groups
  days_to_compare <- alpha_long %>%
    filter(.data[[group_var]] != "untreated") %>%
    distinct(.data[[day_var]]) %>%
    pull()

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
