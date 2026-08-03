#############################################################
## Prepare annotations for ggplot
#############################################################

prepare_alpha_annotations <- function(alpha_stats){

  annotations <- data.frame()

  for(name in names(alpha_stats$details)){

    res <- alpha_stats$details[[name]]

    #########################################################
    ## ANOVA
    #########################################################

    if(res$test == "One-way ANOVA"){

      posthoc <- as.data.frame(res$posthoc)

      groups <- do.call(
        rbind,
        strsplit(posthoc$contrast, " - ")
      )

      posthoc$group1 <- groups[,1]
      posthoc$group2 <- groups[,2]

      posthoc$p.adj <- posthoc$p.value

      posthoc$p.adj.signif <- dplyr::case_when(
        posthoc$p.adj <= 0.0001 ~ "****",
        posthoc$p.adj <= 0.001  ~ "***",
        posthoc$p.adj <= 0.01   ~ "**",
        posthoc$p.adj <= 0.05   ~ "*",
        TRUE ~ "ns"
      )

    } else {

      #######################################################
      ## Kruskal
      #######################################################

      posthoc <- res$posthoc

    }

    #########################################################
    ## Add metadata
    #########################################################

    posthoc$Index <- res$index
    posthoc$Day <- res$day

    annotations <- dplyr::bind_rows(
      annotations,
      posthoc
    )

  }

  return(annotations)

}
