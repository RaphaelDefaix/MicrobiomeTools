#############################################################
## 1. Header
#############################################################
#
#  MicrobiomeTools
#
#  Script      : Composition_plot.R
#  Author      : Raphael Defaix
#  Repository  : https://github.com/RaphaelDefaix/MicrobiomeTools
#
#  Version     : 0.1.0
#  Created     : July 2026
#  Description :
#
#  Create publication-quality stacked barplots of microbial
#  composition from phyloseq objects.
#
#############################################################
## CHANGELOG
#############################################################

# v0.1.0
#
# - First reference script
# - Manual colour palette
# - Top taxa selection
# - Publication-quality figure

#############################################################
## Required input
#############################################################

## This script expects:
##
## - A cleaned phyloseq object named 'ps1'
##
## The object should already have:
##
## • chloroplasts removed
## • mitochondria removed
## • contaminants removed (if applicable)
##

#############################################################
## 2. Packages
#############################################################

library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library(scales)
library(stringr)
library(tidyr)

#############################################################
## Modify only the parameters below
#############################################################
#############################################################
## 3. Parameters
#############################################################

## Taxonomic level
  tax_level <- "Genus"
  
## Number of most abundant taxa to display.
## Remaining taxa are merged into "Other".
top_taxa <- 15
  
## Sample metadata
  group_var <- "group"
  day_var <- "day"

#############################################################
#############################################################

## Select one experiment only
## Example:
## experiment <- "RD09"
##
## Use NULL to analyse all samples

experiment <- NULL

## Export options
figure_width  <- 8
figure_height <- 5
figure_dpi    <- 600

export_figure <- TRUE

## Figure title
figure_title <- "Genus composition"

## Output filename
output_name <- "Composition_plot"

output_directory <- "Figures"

#############################################################
## End of user parameters
#############################################################

#############################################################
## 4. Data preparation
#############################################################

## Select experiment (if specified)

if (!is.null(experiment)) {

  ps <- subset_samples(
    ps1,
    exp == experiment
  )

} else {

  ps <- ps1

}

## Remove taxa with zero counts

ps <- prune_taxa(
  taxa_sums(ps) > 0,
  ps
)

## Transform counts to relative abundance

ps.rel <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)

## Aggregate taxa

ps.rel <- tax_glom(
  ps.rel,
  taxrank = tax_level,
  NArm = FALSE
)

## Convert phyloseq object to data frame

df <- psmelt(ps.rel)

#############################################################
## Create Taxon column
#############################################################

df$Taxon <- case_when(

  !is.na(df[[tax_level]]) ~ df[[tax_level]],

  !is.na(df$Family) ~ paste0("Unclassified_", df$Family),

  !is.na(df$Order) ~ paste0("Unclassified_", df$Order),

  !is.na(df$Class) ~ paste0("Unclassified_", df$Class),

  !is.na(df$Phylum) ~ paste0("Unclassified_", df$Phylum),
  
  TRUE ~ "Unknown"

)

## Check required columns

required_columns <- c(
  "Sample",
  "Abundance",
  group_var,
  day_var,
  tax_level
)

missing_columns <- setdiff(
  required_columns,
  colnames(df)
)

if (length(missing_columns) > 0) {

  stop(
    paste(
      "Missing columns:",
      paste(missing_columns,
            collapse = ", ")
    )
  )

}

#############################################################
## 5. Mean abundance calculation
#############################################################

## Sum abundance within each sample

df_sample <- df %>%
  group_by(
    Sample,
    .data[[day_var]],
    .data[[group_var]],
    Taxon
  ) %>%
  summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  )

## Calculate mean abundance for each day and group

df_mean <- df_sample %>%
  group_by(
    .data[[day_var]],
    .data[[group_var]],
    Taxon
  ) %>%
  summarise(
    MeanAbundance = mean(Abundance),
    .groups = "drop"
  )


#############################################################
## 6. Top taxa selection
#############################################################

## Calculate the total abundance of each taxon

top_taxa_list <- df_mean %>%
  group_by(Taxon) %>%
  summarise(
    TotalAbundance = sum(MeanAbundance),
    .groups = "drop"
  ) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = top_taxa) %>%
  pull(Taxon)

## Merge remaining taxa into "Other"

df_plot <- df_mean %>%
  mutate(
    Taxon = if_else(
      Taxon %in% top_taxa_list,
      Taxon,
      "Other"
    )
  ) %>%
  group_by(
    .data[[day_var]],
    .data[[group_var]],
    Taxon
  ) %>%
  summarise(
    MeanAbundance = sum(MeanAbundance),
    .groups = "drop"
  )

## Calculate the final abundance table

taxa_order <- df_plot %>%
  group_by(Taxon) %>%
  summarise(
    TotalAbundance = sum(MeanAbundance),
    .groups = "drop"
  ) %>%
  arrange(desc(TotalAbundance)) %>%
  pull(Taxon)

## Always display "Other" last

taxa_order <- c(
  setdiff(taxa_order, "Other"),
  "Other"
)

## Convert Taxon to ordered factor

df_plot$Taxon <- factor(
  df_plot$Taxon,
  levels = rev(taxa_order)
)

#############################################################
## 7. Colour palette
#############################################################

## Default colour palette

taxon_colors <- c(

  "Unclassified_Muribaculaceae"     = "#D55E00",
  "Muribaculum"                     = "#E69F00",
  "Parasutterella"                  = "#00BFC4",
  "Bacteroides"                     = "#0072B2",
  "Parabacteroides"                 = "#56B4E9",
  "Prevotellaceae UCG-001"          = "#009E73",
  "Alistipes"                       = "#66A61E",

  "Akkermansia"                     = "#CC79A7",

  "Blautia"                         = "#7570B3",
  "Lachnospiraceae NK4A136 group"   = "#8DA0CB",
  "Colidextribacter"                = "#A6761D",

  "Escherichia-Shigella"            = "#E41A1C",
  "Unclassified_Enterobacteriaceae" = "#FB8072",
  "Proteus"                         = "#A50F15",

  "Dubosiella"                      = "#1B9E77",

  "Other"                           = "grey80"

)

## Keep only colours corresponding to displayed taxa

taxon_colors <- taxon_colors[
  levels(df_plot$Taxon)
]

#############################################################
## 8. Create figure
#############################################################

composition_plot <- ggplot(
  df_plot,
  aes(
    x = .data[[day_var]],
    y = MeanAbundance,
    fill = Taxon
  )
) +

  geom_col(
    colour = "black",
    linewidth = 0.2
  ) +

  facet_wrap(
    vars(.data[[group_var]]),
    nrow = 1
  ) +

  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +

  scale_fill_manual(
    values = taxon_colors
  ) +

  labs(
    title = figure_title,
    x = "Day",
    y = "Relative abundance",
    fill = tax_level
  ) +

  theme_bw() +

  theme(

    panel.grid = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),

    strip.background = element_rect(
      fill = "grey90",
      colour = "black"
    ),

    strip.text = element_text(
      face = "bold"
    ),

    legend.title = element_text(
      face = "bold"
    ),

    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )

  )

#############################################################
## 9. Export figure
#############################################################

if (export_figure) {

  ## Create output directory if it does not exist

  dir.create(
    output_directory,
    showWarnings = FALSE,
    recursive = TRUE
  )

  ## Export figure

  ggsave(
    filename = file.path(
      output_directory,
      paste0(output_name, ".png")
    ),
    plot = composition_plot,
    width = figure_width,
    height = figure_height,
    dpi = figure_dpi
  )

}
