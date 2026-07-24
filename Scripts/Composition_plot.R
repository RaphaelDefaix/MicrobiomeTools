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
## 2. User parameters
#############################################################

## Taxonomic level
tax_level <- "Genus"

## Number of taxa displayed
top_taxa <- 15

## Sample metadata
group_var <- "group"
day_var <- "day"

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
## 3. Data preparation
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

## Remove empty taxa

ps <- prune_taxa(taxa_sums(ps) > 0, ps)

## Transform counts to relative abundance

ps.rel <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)

#############################################################
## 4. Mean abundance calculation
#############################################################

## Calculate mean abundance per group

#############################################################
## 5. Top taxa selection
#############################################################

## Select Top taxa

## Merge remaining taxa into "Other"

#############################################################
## 6. Colour palette
#############################################################

## Manual palette (v0.1)

#############################################################
## 7. Create figure
#############################################################

## ggplot

#############################################################
## 8. Export figure
#############################################################

## Save figure
