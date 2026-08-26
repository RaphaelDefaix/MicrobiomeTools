############################################################
## GENUS ABUNDANCE — RD03 + RD09 + RD10
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
## Taxonomy
############################################################

tax_df <- as.data.frame(
    tax_table(ps.rel),
    stringsAsFactors = FALSE
)

############################################################
## Create taxon names
############################################################

tax_df$Taxon <- ifelse(
    !is.na(tax_df$Genus) & tax_df$Genus != "",
    as.character(tax_df$Genus),
    paste0(
        "Unclassified_",
        as.character(tax_df$Family)
    )
)

tax_df$Taxon <- gsub(
    " ",
    "_",
    tax_df$Taxon
)

tax_df$Taxon <- make.unique(
    tax_df$Taxon
)

rownames(tax_df) <- taxa_names(ps.rel)

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

otu_df$Taxon <- tax_df[
    rownames(otu_df),
    "Taxon"
]

############################################################
## Aggregate taxa
############################################################

taxon_abundance <- otu_df %>%
    group_by(Taxon) %>%
    summarise(
        across(
            everything(),
            sum,
            na.rm = TRUE
        )
    )

############################################################
## Convert to long format
############################################################

df_long <- taxon_abundance %>%
    pivot_longer(
        cols = -Taxon,
        names_to = "SampleID",
        values_to = "Abundance"
    )

############################################################
## Add metadata
############################################################

metadata <- data.frame(
    sample_data(ps.rel),
    stringsAsFactors = FALSE
)

metadata$SampleID <- rownames(metadata)

df_long <- df_long %>%
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

df_long <- df_long %>%
    filter(
        exp %in% c(
            "RD03",
            "RD09",
            "RD10"
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
    "Lachnospiraceae_NK4A136_group",
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

