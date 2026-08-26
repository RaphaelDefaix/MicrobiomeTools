############################################################
## TOP 6 PHYLA — RD03 + RD09 + RD10
############################################################

library(phyloseq)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)

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
## Taxonomy — Phylum
############################################################

tax_df <- as.data.frame(
    tax_table(ps.rel),
    stringsAsFactors = FALSE
)

tax_df$Phylum <- as.character(tax_df$Phylum)

############################################################
## Abundance matrix
############################################################

otu_df <- as.data.frame(
    otu_table(ps.rel)
)

if (!taxa_are_rows(ps.rel)) {
    otu_df <- t(otu_df)
}

otu_df <- as.data.frame(otu_df)

otu_df$Phylum <- tax_df[
    rownames(otu_df),
    "Phylum"
]

############################################################
## Remove missing Phylum
############################################################

otu_df <- otu_df %>%
    filter(
        !is.na(Phylum),
        Phylum != ""
    )

############################################################
## Aggregate ASVs by Phylum
############################################################

phylum_abundance <- otu_df %>%
    group_by(Phylum) %>%
    summarise(
        across(
            everything(),
            sum,
            na.rm = TRUE
        )
    )

############################################################
## Long format
############################################################

df_phylum <- phylum_abundance %>%
    pivot_longer(
        cols = -Phylum,
        names_to = "SampleID",
        values_to = "Abundance"
    )

############################################################
## Metadata
############################################################

metadata <- data.frame(
    sample_data(ps.rel),
    stringsAsFactors = FALSE
)

metadata$SampleID <- rownames(metadata)

df_phylum <- df_phylum %>%
    left_join(
        metadata %>%
            select(
                SampleID,
                exp,
                day,
                group
            ),
        by = "SampleID"
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

top6_phyla <- df_phylum %>%
    group_by(Phylum) %>%
    summarise(
        MeanAbundance = mean(
            Abundance,
            na.rm = TRUE
        )
    ) %>%
    arrange(
        desc(MeanAbundance)
    ) %>%
    slice_head(
        n = 6
    )

top6_phyla

############################################################
## TOP 4 PHYLA
############################################################

top4_phyla <- df_phylum %>%
    group_by(Phylum) %>%
    summarise(
        MeanAbundance = mean(
            Abundance,
            na.rm = TRUE
        )
    ) %>%
    arrange(
        desc(MeanAbundance)
    ) %>%
    slice_head(n = 4)

top4_phyla

############################################################
## Keep TOP 4
############################################################

df_top4_phyla <- df_phylum %>%
    filter(
        Phylum %in% top4_phyla$Phylum
    )

############################################################
## Update Phylum names
############################################################

df_top4_phyla <- df_top4_phyla %>%
    mutate(
        Phylum = recode(
            Phylum,
            "Firmicutes" = "Bacillota",
            "Proteobacteria" = "Pseudomonadota"
        )
    )

df_top4_phyla$Phylum <- factor(
    df_top4_phyla$Phylum,
    levels = c(
        "Bacteroidota",
        "Pseudomonadota",
        "Verrucomicrobiota",
        "Bacillota"
    )
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
