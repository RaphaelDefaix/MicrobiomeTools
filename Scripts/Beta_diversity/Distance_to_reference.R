
#############################################################
## Libraries
#############################################################

library(phyloseq)
library(dplyr)
library(vegan)
library(emmeans)
library(car)
library(FSA)
library(ggplot2)


#############################################################
## Input phyloseq object
#############################################################

ps <- ps_RD03050910

#############################################################
## Distance
#############################################################

distance_method <- "bray"

#############################################################
## Reference group
#############################################################

reference_group <- "Untreated"

#############################################################
## Variable containing the groups
#############################################################

group_variable <- "group"

#############################################################
## Optional filtering
#############################################################

if (!inherits(ps, "phyloseq")) {
    stop("'ps' must be a phyloseq object.")
}

if (!group_variable %in% colnames(sample_data(ps))) {
    stop("'group_variable' is not present in sample_data(ps).")
}

#############################################################
## Check UniFrac requirements
#############################################################

if (distance_method %in% c("unifrac", "wunifrac")) {

    if (is.null(phy_tree(ps, errorIfNULL = FALSE))) {
        stop("A phylogenetic tree is required for UniFrac distances.")
    }

}

#############################################################
## Optional filtering
#############################################################

experiment <- NULL

exclude_cages <- NULL

exclude_groups <- NULL

exclude_days <- NULL

#############################################################
## Optional condition filtering
#############################################################

exclude_conditions <- NULL

## Examples:
## exclude_conditions <- data.frame(
##     group = c("Antibiotic","Antibiotic","Antibiotic"),
##     day   = c("d1","d7","d13")
## )

#############################################################
## Optional experiment filtering
#############################################################

if (!is.null(experiment)) {

    ps <- subset_samples(
        ps,
        exp %in% experiment
    )

}

#############################################################
## Optional cage filtering
#############################################################

if (!is.null(exclude_cages)) {

    ps <- subset_samples(
        ps,
        !(cage %in% exclude_cages)
    )

}

#############################################################
## Optional group filtering
#############################################################

if (!is.null(exclude_groups)) {

    ps <- subset_samples(
        ps,
        !(group %in% exclude_groups)
    )

}

#############################################################
## Optional day filtering
#############################################################

if (!is.null(exclude_days)) {

    ps <- subset_samples(
        ps,
        !(day %in% exclude_days)
    )

}

#############################################################
## Optional condition filtering
#############################################################

if (!is.null(exclude_conditions)) {

    metadata <- data.frame(sample_data(ps))

    keep <- rep(TRUE, nrow(metadata))

    for (i in seq_len(nrow(exclude_conditions))) {

        remove <- rep(TRUE, nrow(metadata))

        for (variable in names(exclude_conditions)) {

            remove <- remove &
                metadata[[variable]] == exclude_conditions[i, variable]

        }

        keep <- keep & !remove

    }

    ps <- prune_samples(keep, ps)

}

#############################################################
## Distance matrix
#############################################################

dist_matrix <- phyloseq::distance(
    ps,
    method = distance_method
)


#############################################################
## Reference samples
#############################################################

metadata <- data.frame(sample_data(ps))

reference_samples <- rownames(
    metadata[
        metadata[[group_variable]] == reference_group,
        ,
        drop = FALSE
    ]
)

if (length(reference_samples) == 0) {
    stop("No samples found for the selected reference group.")
}


#############################################################
## Distance to reference
#############################################################

dist_df <- as.matrix(dist_matrix)

distance_table <- data.frame(
    SampleID = rownames(dist_df),
    Distance_to_reference = NA_real_
)

for (i in seq_len(nrow(distance_table))) {

    sample_id <- distance_table$SampleID[i]

    ## Reference samples to use
    refs <- reference_samples

    ## If the sample belongs to the reference group,
    ## remove it from the reference set
    if (sample_id %in% reference_samples) {
        refs <- setdiff(reference_samples, sample_id)
    }

    distance_table$Distance_to_reference[i] <-
        mean(
            dist_df[sample_id, refs],
            na.rm = TRUE
        )

}


#############################################################
## Add metadata
#############################################################

metadata <- data.frame(sample_data(ps))

metadata$SampleID <- rownames(metadata)

distance_table <- dplyr::left_join(
    distance_table,
    metadata,
    by = "SampleID"
)

#############################################################
## Order groups
#############################################################

group_order <- c(
    "Untreated",
    "Antibiotic",
    "WT",
    "lsrK",
    "lsrR"
)

distance_table$group <- factor(
    distance_table$group,
    levels = group_order
)

#############################################################
## Remove reference group for plotting
#############################################################

distance_plot <- dplyr::filter(
    distance_table,
    group != reference_group
)

#############################################################
## Statistics
#############################################################

statistics <- "anova"

## Reserved for future versions:
## "anova", "kruskal", "auto"

#############################################################
## One-way ANOVA
#############################################################

anova_model <- aov(
    Distance_to_reference ~ group,
    data = distance_plot
)

#############################################################
## ANOVA assumptions
#############################################################

shapiro_result <- shapiro.test(
    residuals(anova_model)
)

print(shapiro_result)

levene_result <- leveneTest(
    Distance_to_reference ~ group,
    data = distance_plot
)

print(levene_result)

#############################################################
## ANOVA results
#############################################################

anova_results <- summary(anova_model)

print(anova_results)

#############################################################
## Tukey post-hoc
#############################################################

pairwise_results <- emmeans(
    anova_model,
    pairwise ~ group,
    adjust = "tukey"
)

pairwise_table <- as.data.frame(
    pairwise_results$contrasts
)

print(pairwise_table)


#############################################################
## Plot
#############################################################

ggplot(
    distance_plot,
    aes(
        x = group,
        y = Distance_to_reference,
        fill = group
    )
) +
    geom_boxplot(
        alpha = 0.6,
        colour = "black"
    ) +
    geom_jitter(
        width = 0.15,
        size = 2,
        shape = 21,
        colour = "black"
    ) +
    scale_fill_manual(
    values = c(
        "Untreated" = "grey70",
        "Antibiotic" = "grey40",
        "WT" = "#03c2fe",
        "lsrK" = "#02d62e",
        "lsrR" = "#ff1130"
    )
)+
    labs(
        x = "Group",
        y = "Bray-Curtis distance to reference"
    ) +
    microbiome_theme()

#############################################################
## Export
#############################################################

write.csv(
    distance_table,
    "Distance_to_reference.csv",
    row.names = FALSE
)

#############################################################
## Export statistics
#############################################################

write.csv(
    pairwise_table,
    "Distance_to_reference_statistics.csv",
    row.names = FALSE
)

