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
#  Version     : 1.0.0
#  Created     : July 2026
#  Description :
#
#  Create publication-quality stacked barplots of microbial
#  composition from phyloseq objects.
#
#############################################################
## CHANGELOG
#############################################################

# v1.0.0
#
# - Added microbiome_palette()
# - Added microbiome_theme()
# - Added biological taxon ordering
# - Added parameter validation
# - Added input validation
# - Improved code readability
# - Automatic figure display

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
library(grid)
library(phyloseq)
library(scales)
library(tidyr)

#############################################################
## Load MicrobiomeTools functions
#############################################################

source("R/microbiome_palette.R")
source("R/microbiome_theme.R")

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

## Day order
day_order <- c(
  "d-4",
  "d-1",
  "d1",
  "d7",
  "d13"
)

## Group order
group_order <- c(
  "untreated",
  "WT",
  "lsrK",
  "lsrR"
)

#############################################################
#############################################################

## Select one experiment only
## Example:
## experiment <- "RD09"
##
## Use NULL to analyse all samples

experiment <- NULL

## Export options
figure_width  <- 12
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
## Validate user parameters
#############################################################

valid_tax_levels <- c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

if (!tax_level %in% valid_tax_levels) {

  stop(
    paste(
      "Invalid tax_level:",
      tax_level
    )
  )

}

#############################################################
## Check required input
#############################################################

if (!exists("ps1")) {

  stop(
    "The phyloseq object 'ps1' was not found."
  )

}

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

if (nsamples(ps) == 0) {

  stop(
    "No samples found after filtering."
  )

}

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

df <- df %>%
  mutate(

    Taxon = case_when(

      !is.na(.data[[tax_level]]) ~ .data[[tax_level]],

      !is.na(Family) ~ paste0("Unclassified_", Family),

      !is.na(Order) ~ paste0("Unclassified_", Order),

      !is.na(Class) ~ paste0("Unclassified_", Class),

      !is.na(Phylum) ~ paste0("Unclassified_", Phylum),

      TRUE ~ "Unknown"

    )

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

#############################################################
## Taxon - Phylum correspondence
#############################################################

taxon_phylum <- tibble(

  Taxon = c(

    "Unclassified_Muribaculaceae",
    "Muribaculum",
    "Bacteroides",
    "Parabacteroides",
    "Prevotellaceae UCG-001",
    "Alistipes",

    "Akkermansia",

    "Lachnospiraceae NK4A136 group",
    "Blautia",
    "Dubosiella",

    "Escherichia-Shigella",
    "Proteus",
    "Unclassified_Enterobacteriaceae"

  ),

  Phylum = c(

    "Bacteroidota",
    "Bacteroidota",
    "Bacteroidota",
    "Bacteroidota",
    "Bacteroidota",
    "Bacteroidota",

    "Verrucomicrobiota",

    "Bacillota",
    "Bacillota",
    "Bacillota",

    "Pseudomonadota",
    "Pseudomonadota",
    "Pseudomonadota"

  )

)

#############################################################
## Phylum order
#############################################################

phylum_order <- c(

  "Bacteroidota",

  "Bacillota",

  "Pseudomonadota",

  "Verrucomicrobiota"

)

#############################################################
## Biological ordering of taxa
#############################################################

taxa_order <- df_plot %>%

  left_join(
    taxon_phylum,
    by = "Taxon"
  ) %>%

  mutate(
    Phylum = factor(
      Phylum,
      levels = phylum_order
    )
  ) %>%

  group_by(
    Phylum,
    Taxon
  ) %>%

  summarise(
    TotalAbundance = sum(MeanAbundance),
    .groups = "drop"
  ) %>%

  arrange(
    Phylum,
    desc(TotalAbundance)
  ) %>%

  pull(Taxon)

taxa_order <- c(
  setdiff(taxa_order, "Other"),
  "Other"
)

df_plot$Taxon <- factor(
  df_plot$Taxon,
  levels = rev(taxa_order)
)

#############################################################
## 7. Colour palette
#############################################################
#############################################################
## Official MicrobiomeTools palette
#############################################################

taxon_colors <- microbiome_palette()

## Keep only colours corresponding to displayed taxa

taxon_colors <- taxon_colors[
  levels(df_plot$Taxon)
]

#############################################################
## 8. Create figure
#############################################################

df_plot[[day_var]] <- factor(
  df_plot[[day_var]],
  levels = day_order
)

df_plot[[group_var]] <- factor(
  df_plot[[group_var]],
  levels = group_order
)
  
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
  linewidth = 0.2,
  width = 0.9
) +

facet_wrap(
  vars(.data[[group_var]]),
  nrow = 1,
  scales = "free_x"
) +

scale_x_discrete(drop = TRUE) +

scale_y_continuous(
  labels = scales::percent_format(accuracy = 1),
  expand = expansion(mult = c(0, 0.02))
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

microbiome_theme()
  
  composition_plot

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
