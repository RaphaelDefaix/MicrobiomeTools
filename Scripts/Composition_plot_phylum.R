#############################################################
## Input object
#############################################################

ps3 <- physeq_RD10

#############################################################

## Select one experiment only
## Example:
## experiment <- "RD10"
##
## Use NULL to analyse all samples

experiment <- NULL

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
  tax_level <- "Phylum"
  
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

## Export options
figure_width  <- 12
figure_height <- 5
figure_dpi    <- 600

export_figure <- TRUE

## Figure title
figure_title <- paste(
  "Phylum composition",
  ifelse(is.null(experiment), "", experiment)
)

## Output filename
output_name <- paste0(
  "Phylum_composition_",
  ifelse(is.null(experiment),"All",experiment)
)

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
    ps3,
    exp == experiment
  )

} else {

  ps <- ps3

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
    Taxon = ifelse(
      is.na(.data[[tax_level]]),
      "Unclassified",
      as.character(.data[[tax_level]])
    )
  ) %>%
  mutate(
    Taxon = recode(
      Taxon,
      "Firmicutes"       = "Bacillota",
      "Proteobacteria"   = "Pseudomonadota",
      "Actinobacteriota" = "Actinomycetota", 
      "Bacteroidetes"      = "Bacteroidota",
      "Epsilonbacteraeota" = "Campylobacterota"
    )
  )
#############################################################
## Create Taxon column
#############################################################

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
## Biological order of phyla
#############################################################

phylum_order <- c(

"Bacteroidota",

"Bacillota",

"Pseudomonadota",

"Verrucomicrobiota",

"Actinomycetota",

"Desulfobacterota",

"Campylobacterota",

"Fusobacteriota",

"Patescibacteria",

"Deferribacterota",

"Unclassified",

"Other"

)

df_plot$Taxon <- factor(
  df_plot$Taxon,
  levels = rev(phylum_order)
)
#############################################################
## 7. Colour palette
#############################################################
#############################################################
## Official MicrobiomeTools palette
#############################################################

phylum_palette <- c(

"Bacteroidota"      ="#2171B5",

"Bacillota"         ="#31A354",

"Pseudomonadota"    ="#CB181D",

"Verrucomicrobiota" ="#F16913",

"Actinomycetota"    ="#756BB1",

"Desulfobacterota"  ="#8C6D31",

"Deferribacterota"  ="#E67E22",

"Unclassified"      ="grey70",

"Other"             ="grey90"

)

## Keep only colours corresponding to displayed taxa

taxon_colors <- phylum_palette[
    levels(df_plot$Taxon)
]

if(any(is.na(taxon_colors))){

  warning(
    "Some displayed phyla have no defined colour."
  )

  print(levels(df_plot$Taxon)[is.na(taxon_colors)])

}

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
