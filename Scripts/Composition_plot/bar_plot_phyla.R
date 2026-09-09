############################################################
## TOP PHYLA — RD03 + RD09 + RD10
############################################################

library(phyloseq)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)

############################################################
## MicrobiomeTools functions
############################################################

source("R/taxonomy_database.R")

############################################################
## Input
############################################################

ps <- ps_RD030910

############################################################
## Relative abundance
############################################################

ps.rel <- transform_sample_counts(
    ps,
    function(x) x / sum(x)
)

############################################################
## Aggregate taxa at Phylum level
############################################################

ps.rel <- tax_glom(
    ps.rel,
    taxrank = "Phylum",
    NArm = FALSE
)

############################################################
## Convert phyloseq to data frame
############################################################

df_phylum <- psmelt(
    ps.rel
)

############################################################
## Standardize Phylum names
############################################################

df_phylum$Phylum <- standardize_phylum(
    df_phylum$Phylum
)


############################################################
## Keep d1, d7 and d13 only
############################################################

df_phylum <- df_phylum %>%
    filter(
        day %in% c(
            "d1",
            "d7",
            "d13"
        )
    )

############################################################
## Keep WT, lsrK and lsrR
############################################################

df_phylum <- df_phylum %>%
    filter(
        group %in% c(
            "WT",
            "lsrK",
            "lsrR"
        )
    )

############################################################
## Order groups
############################################################

df_phylum$group <- factor(
    df_phylum$group,
    levels = c(
        "WT",
        "lsrK",
        "lsrR"
    )
)

############################################################
## Check
############################################################

table(
    df_phylum$exp,
    df_phylum$day,
    df_phylum$group
)

############################################################
## Keep RD03 + RD09 + RD10
############################################################

df_phylum <- df_phylum %>%
    filter(
        exp %in% c(
            "RD03",
            "RD09",
            "RD10"
        )
    ) %>%
    filter(
        !is.na(group)
    )

############################################################
## Calculate mean abundance for each Phylum
############################################################

phylum_abundance <- df_phylum %>%
    group_by(Phylum) %>%
    summarise(
        MeanAbundance = mean(
            Abundance,
            na.rm = TRUE
        ),
        .groups = "drop"
    ) %>%
    arrange(
        desc(MeanAbundance)
    )

############################################################
## TOP 6 PHYLA
############################################################

top6_phyla <- phylum_abundance %>%
    slice_head(
        n = 6
    )

top6_phyla

############################################################
## TOP 4 PHYLA
############################################################

top4_phyla <- phylum_abundance %>%
    slice_head(
        n = 4
    )

top4_phyla

############################################################
## Keep TOP 4
############################################################

df_top4_phyla <- df_phylum %>%
    filter(
        Phylum %in% top4_phyla$Phylum
    )

############################################################
## Order Phyla
############################################################

phylum_order <- c(
    "Bacteroidota",
    "Bacillota",
    "Pseudomonadota",
    "Verrucomicrobiota"
)

df_top4_phyla$Phylum <- factor(
    df_top4_phyla$Phylum,
    levels = phylum_order
)

############################################################
## PLOT — TOP 4 PHYLA
############################################################

plot_top4_phyla <- ggplot(
    df_top4_phyla,
    aes(
        x = group,
        y = Abundance
    )
) +

    geom_boxplot(
        width = 0.55,
        outlier.shape = NA
    ) +

    geom_jitter(
        width = 0.10,
        size = 2,
        alpha = 0.8
    ) +

    facet_wrap(
        ~ Phylum,
        ncol = 2,
        scales = "free_y"
    ) +

    scale_x_discrete(
        limits = c(
            "WT",
            "lsrK",
            "lsrR"
        )
    ) +

    labs(
        title = "Phylum abundance",
        x = NULL,
        y = "Relative abundance"
    ) +

    theme_classic(
        base_size = 14
    ) +

    theme(

        strip.background = element_rect(
            colour = "black",
            fill = "white"
        ),

        strip.text = element_text(
            face = "italic",
            margin = margin(
                t = 5,
                b = 12
            )
        ),

        strip.placement = "outside",

        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        ),

        plot.title = element_text(
            hjust = 0
        ),

        panel.spacing = unit(
            1.2,
            "lines"
        )
    ) +

    scale_y_continuous(
        expand = expansion(
            mult = c(0.05, 0.25)
        )
    )

plot_top4_phyla

############################################################
## ONE-WAY ANOVA + TUKEY — TOP 4 PHYLA
############################################################

anova_phylum <- function(df, phylum_name) {
    
    data_phylum <- df %>%
        filter(
            Phylum == phylum_name,
            !is.na(Abundance),
            !is.na(group)
        )
    
    data_phylum$group <- factor(
        data_phylum$group,
        levels = c(
            "WT",
            "lsrK",
            "lsrR"
        )
    )
    
    ########################################################
    ## One-way ANOVA
    ########################################################
    
    model <- aov(
        Abundance ~ group,
        data = data_phylum
    )
    
    ########################################################
    ## Tukey HSD
    ########################################################
    
    tukey <- TukeyHSD(model)
    
    ########################################################
    ## Print
    ########################################################
    
    cat("\n========================================\n")
    cat("Phylum:", phylum_name, "\n")
    cat("========================================\n\n")
    
    cat("One-way ANOVA:\n")
    print(summary(model))
    
    cat("\nTukey HSD:\n")
    print(tukey)
    
    ########################################################
    ## Return results
    ########################################################
    
    return(
        list(
            model = model,
            anova = summary(model),
            tukey = tukey
        )
    )
}

############################################################
## Run ANOVA + Tukey for all 4 Phyla
############################################################

results_phyla <- lapply(
    levels(df_top4_phyla$Phylum),
    function(x) {
        anova_phylum(
            df_top4_phyla,
            x
        )
    }
)

names(results_phyla) <- levels(
    df_top4_phyla$Phylum
)
