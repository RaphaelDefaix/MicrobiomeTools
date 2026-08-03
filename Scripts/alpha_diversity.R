#############################################################
## Input object
#############################################################

ps3 <- physeq_RD09

## Select one experiment only
## Example:
## experiment <- "RD09"
##
## Use NULL to analyse all samples

experiment <- "RD09"

#############################################################
## Figure title
#############################################################

figure_title <- paste(
  "Alpha diversity",
  ifelse(is.null(experiment), "", experiment)
)

#############################################################
## 1. Header
#############################################################
#
#  MicrobiomeTools
#
#  Script      : Alpha_diversity.R
#  Author      : Raphael Defaix
#  Repository  : https://github.com/RaphaelDefaix/MicrobiomeTools
#
#  Version     : 1.0.0
#
#############################################################

library(phyloseq)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(tidyr)

source("R/export_alpha_statistics.R")
source("R/microbiome_palette.R")
source("R/microbiome_theme.R")
source("R/run_alpha_statistics.R")
source("R/prepare_alpha_annotations.R")

#############################################################
## Parameters
#############################################################

group_var <- "group"

day_var <- "day"

group_order <- c(
  "Untreated",
  "Antibiotic",
  "WT",
  "lsrK",
  "lsrR"
)

day_order <- c(
  "d-4",
  "d-1",
  "d1",
  "d7",
  "d13"
)

#############################################################
## Check input
#############################################################

if (!inherits(ps3, "phyloseq")) {

  stop(
    "'ps3' must be a phyloseq object."
  )

}

#############################################################
## Select experiment
#############################################################

if (!is.null(experiment)) {

  ps <- subset_samples(
    ps3,
    exp == experiment
  )

} else {

  ps <- ps3

}

#############################################################
## Alpha diversity indices
#############################################################

alpha_indices <- c(
  "Observed",
  "Shannon",
  "InvSimpson"
)


## Note:
## estimate_richness() may warn if singletons have already
## been removed. This mainly affects Chao1 estimates.

#############################################################
## Calculate alpha diversity
#############################################################

alpha_df <- estimate_richness(
  ps,
  measures = alpha_indices
)

############################################################
## Add metadata
#############################################################

metadata <- data.frame(
  sample_data(ps)
)

metadata$Sample <- rownames(metadata)

alpha_df$Sample <- rownames(alpha_df)

alpha_df <- left_join(
  alpha_df,
  metadata,
  by = "Sample"
)

alpha_df[[group_var]] <- factor(
  alpha_df[[group_var]],
  levels = group_order
)

alpha_df[[day_var]] <- factor(
  alpha_df[[day_var]],
  levels = day_order
)

#############################################################
## Convert to long format
#############################################################

alpha_long <- alpha_df %>%

  pivot_longer(

    cols = all_of(alpha_indices),

    names_to = "Index",

    values_to = "Value"

  )


#############################################################
## Statistical analysis
#############################################################

alpha_stats <- run_alpha_statistics(
    alpha_long = alpha_long,
    group_var = group_var,
    day_var = day_var,
    alpha_indices = alpha_indices
)

#############################################################
## Display statistical summary
#############################################################

cat("\n")
cat("=========================================\n")
cat(" Alpha diversity statistical summary\n")
cat("=========================================\n\n")

print(alpha_stats$summary)

#############################################################
## Export statistical results
#############################################################
#############################################################
## Prepare annotations
#############################################################

annotations <- prepare_alpha_annotations(alpha_stats)

annotations$comparisons <- Map(
  c,
  annotations$group1,
  annotations$group2
)

#############################################################
## Export statistics
#############################################################

annotations_export <- annotations %>%
  dplyr::select(-comparisons)

export_alpha_statistics(
    summary = alpha_stats$summary,
    posthoc = annotations_export,
    experiment = experiment
)

#############################################################
## Order indices
#############################################################

alpha_long$Index <- factor(

  alpha_long$Index,

  levels = alpha_indices

)


#############################################################
## Check table
#############################################################

head(alpha_df)

str(alpha_df)

summary(alpha_df)


#############################################################
## Plot alpha diversity
#############################################################

group_colors <- microbiome_group_palette()

alpha_plot <- ggplot(

  alpha_long,

  aes(

    x = .data[[day_var]],

    y = Value,

    fill = .data[[group_var]]

  )

) +

geom_boxplot(

  aes(group = interaction(.data[[day_var]], .data[[group_var]])),

  alpha = 0.6,

  colour = "black",

  outlier.shape = NA,

  position = position_dodge(width = 0.75),

  width = 0.65

) +

geom_jitter(

  aes(

    colour = .data[[group_var]],

    group = .data[[group_var]]

  ),

  position = position_jitterdodge(

    jitter.width = 0.15,

    dodge.width = 0.75

  ),

  size = 1.8,

  alpha = 0.9

) +

facet_wrap(
  ~Index,
  nrow = 1,
  scales = "free_y"
  
) +


scale_fill_manual(

  values = group_colors

) +

scale_colour_manual(

  values = group_colors

) +

labs(

  title = figure_title,

  x = "Day",

  y = NULL

)+

microbiome_theme()

alpha_plot <- add_alpha_annotations(
    plot = alpha_plot,
    annotations = annotations,
    alpha_long = alpha_long
)

alpha_plot
