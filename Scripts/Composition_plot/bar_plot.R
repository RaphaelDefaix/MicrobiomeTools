#############################################################
## Libraries
#############################################################

library(ggplot2)
library(dplyr)
library(tidyr)


#############################################################
## Input data
#############################################################

data <- read.csv(
    "RD03050910_Genus_relative_abundance.csv",
    check.names = FALSE
)


#############################################################
## Taxa to plot
#############################################################

taxa_to_plot <- c(
    "Unclassified_Muribaculaceae",
    "Bacteroides",
    "Akkermansia",
    "Parabacteroides",
    "Lachnospiraceae NK4A136 group",
    "Unclassified_Enterobacteriaceae",
    "Blautia",
    "Escherichia-Shigella"
)


#############################################################
## Filter samples
#############################################################

data_plot <- data %>%
    
    filter(
        exp %in% c("RD03", "RD09"),
        day %in% c("d1", "d7", "d13"),
        group %in% c("WT", "lsrK", "lsrR")
    )


#############################################################
## Convert to long format
#############################################################

data_long <- data_plot %>%
    
    pivot_longer(
        cols = all_of(taxa_to_plot),
        names_to = "Taxon",
        values_to = "Abundance"
    )


#############################################################
## Order groups
#############################################################

data_long$group <- factor(
    data_long$group,
    levels = c(
        "WT",
        "lsrK",
        "lsrR"
    )
)


#############################################################
## Boxplots
#############################################################

ggplot(
    data_long,
    aes(
        x = group,
        y = Abundance
    )
) +
    
    geom_boxplot(
        width = 0.6,
        outlier.shape = NA
    ) +
    
    geom_jitter(
        width = 0.12,
        size = 2,
        alpha = 0.7
    ) +
    
    facet_wrap(
        ~ Taxon,
        scales = "free_y",
        ncol = 4
    ) +
    
    theme_classic(
        base_size = 14
    ) +
    
    labs(
        title = "Genus abundance",
        x = NULL,
        y = "Relative abundance"
    ) +
    
    theme(
        strip.text = element_text(
            face = "italic",
            size = 11
        ),
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        )
    )
