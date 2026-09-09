############################################################
## GENUS ABUNDANCE — RD03 + RD09 + RD10
############################################################

library(phyloseq)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)

############################################################
## MicrobiomeTools functions
############################################################

source("R/create_taxon_names.R")

############################################################
## Input
############################################################

ps <- ps_RD030910

############################################################
## Groups
############################################################

sample_data(ps)$group <- factor(
    sample_data(ps)$group,
    levels = c("WT", "lsrK", "lsrR")
)

############################################################
## Relative abundance
############################################################

ps.rel <- transform_sample_counts(
    ps,
    function(x) x / sum(x)
)

############################################################
## Aggregate taxa at Genus level
############################################################

ps.rel <- tax_glom(
    ps.rel,
    taxrank = "Genus",
    NArm = FALSE
)

############################################################
## Convert phyloseq to data frame
############################################################

df_long <- psmelt(
    ps.rel
)

############################################################
## Create standardized taxon names
############################################################

df_long$Taxon <- create_taxon_names(
    df_long,
    "Genus"
)


############################################################
## Keep RD03 + RD09 + RD10
############################################################

df_long <- df_long %>%
    filter(
        exp %in% c(
            "RD03",
            "RD09",
            "RD10"
        )
    )
############################################################
## Keep d1, d7 and d13 only
############################################################

df_long <- df_long %>%
    filter(
        day %in% c(
            "d1",
            "d7",
            "d13"
        )
    )

############################################################
## Remove missing groups
############################################################

df_long <- df_long %>%
    filter(
        !is.na(group)
    )
############################################################
## Check metadata
############################################################

df_long %>%
    select(
        Sample,
        exp,
        day,
        group,
        Taxon,
        Abundance
    ) %>%
    head()

############################################################
## Check
############################################################

table(
    df_long$group
)

taxa_8 <- c(
    "Akkermansia",
    "Bacteroides",
    "Blautia",
    "Escherichia-Shigella",
    "Lachnospiraceae NK4A136 group",
    "Parabacteroides",
    "Unclassified_Enterobacteriaceae",
    "Unclassified_Muribaculaceae"
)

df_8taxa <- df_long %>%
    filter(
        Taxon %in% taxa_8
    )
setdiff(
    taxa_8,
    unique(df_long$Taxon)
)
############################################################
## Plot
############################################################

plot_8taxa <- ggplot(
    df_8taxa,
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
    ~ Taxon,
    ncol = 4,
    scales = "free_y"
) +

scale_y_continuous(
    expand = expansion(
        mult = c(0.02, 0.20)
    )

)+

    scale_x_discrete(
        limits = c(
            "WT",
            "lsrK",
            "lsrR"
        )
    ) +

    labs(
        title = "Genus abundance",
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
        b = 5
    ),
    vjust = 0.2
),panel.spacing.y = unit(0.8, "cm"),

        strip.placement = "outside",

        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        ),

        plot.title = element_text(
            hjust = 0
        )
    )

plot_8taxa


## ANOVA ##

############################################################
## One-way ANOVA function
############################################################

anova_taxon <- function(
    df,
    taxon_name
) {

    data_taxon <- df %>%
        filter(
            Taxon == taxon_name
        ) %>%
        filter(
            !is.na(Abundance),
            !is.na(group)
        )

    data_taxon$group <- factor(
        data_taxon$group,
        levels = c(
            "WT",
            "lsrK",
            "lsrR"
        )
    )

    ########################################################
    ## ANOVA
    ########################################################

    model <- aov(
        Abundance ~ group,
        data = data_taxon
    )

    anova_table <- summary(model)

    ########################################################
    ## Tukey
    ########################################################

    tukey <- TukeyHSD(
        model
    )

    ########################################################
    ## Print
    ########################################################

    cat("\n========================================\n")
    cat("Taxon:", taxon_name, "\n")
    cat("========================================\n\n")

    cat("One-way ANOVA:\n")

    print(
        anova_table
    )

    cat("\nTukey HSD:\n")

    print(
        tukey
    )

    ########################################################
    ## Return results
    ########################################################

    return(
        list(
            model = model,
            anova = anova_table,
            tukey = tukey
        )
    )
}

############################################################
## Run ANOVA for all 8 taxa
############################################################

anova_results <- lapply(
    taxa_8,
    function(taxon) {
        anova_taxon(
            df_8taxa,
            taxon
        )
    }
)

names(anova_results) <- taxa_8

