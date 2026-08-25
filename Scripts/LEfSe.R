#############################################################
## Input phyloseq object
#############################################################

ps <- ps_RD03050910


#############################################################
## Filtering
#############################################################

experiment <- NULL

exclude_cages <- NULL

exclude_groups <- NULL

exclude_days <- NULL

exclude_conditions <- NULL


#############################################################
## LEfSe parameters
#############################################################

group_var <- "group"

block_var <- NULL

tax_level <- "Genus"

lda_cutoff <- 2

kw_cutoff <- 0.05

wilcox_cutoff <- 0.05

output_directory <- "Results"


#############################################################
## Libraries
#############################################################

library(phyloseq)
library(SummarizedExperiment)
library(S4Vectors)
library(lefser)
library(dplyr)


#############################################################
## MicrobiomeTools functions
#############################################################

source("R/create_taxon_names.R")
source("R/phyloseq_to_se.R")


#############################################################
## Check input
#############################################################

if (!inherits(ps, "phyloseq")) {
    
    stop(
        "Input object must be a phyloseq object."
    )
    
}


#############################################################
## Check group variable
#############################################################

if (!group_var %in% colnames(sample_data(ps))) {
    
    stop(
        paste(
            "Group variable not found:",
            group_var
        )
    )
    
}


#############################################################
## Check block variable
#############################################################

if (!is.null(block_var)) {
    
    if (!block_var %in% colnames(sample_data(ps))) {
        
        stop(
            paste(
                "Block variable not found:",
                block_var
            )
        )
        
    }
    
}


#############################################################
## Experiment filtering
#############################################################

if (!is.null(experiment)) {
    
    ps <- subset_samples(
        ps,
        exp %in% experiment
    )
    
}


#############################################################
## Cage filtering
#############################################################

if (!is.null(exclude_cages)) {
    
    ps <- subset_samples(
        ps,
        !(cage %in% exclude_cages)
    )
    
}


#############################################################
## Group filtering
#############################################################

if (!is.null(exclude_groups)) {
    
    ps <- subset_samples(
        ps,
        !(group %in% exclude_groups)
    )
    
}


#############################################################
## Day filtering
#############################################################

if (!is.null(exclude_days)) {
    
    ps <- subset_samples(
        ps,
        !(day %in% exclude_days)
    )
    
}


#############################################################
## Combination filtering
#############################################################

if (!is.null(exclude_conditions)) {
    
    metadata <- data.frame(
        sample_data(ps)
    )
    
    #########################################################
    ## Check variables
    #########################################################
    
    missing_variables <- setdiff(
        names(exclude_conditions),
        colnames(metadata)
    )
    
    if (length(missing_variables) > 0) {
        
        stop(
            paste(
                "Unknown variable(s):",
                paste(
                    missing_variables,
                    collapse = ", "
                )
            )
        )
        
    }
    
    
    #########################################################
    ## Determine samples to keep
    #########################################################
    
    keep <- rep(
        TRUE,
        nrow(metadata)
    )
    
    for (i in seq_len(
        nrow(exclude_conditions)
    )) {
        
        remove <- rep(
            TRUE,
            nrow(metadata)
        )
        
        for (variable in names(
            exclude_conditions
        )) {
            
            remove <- remove &
                metadata[[variable]] ==
                exclude_conditions[i, variable]
            
        }
        
        keep <- keep & !remove
        
    }
    
    ps <- prune_samples(
        keep,
        ps
    )
    
}


#############################################################
## Remove empty samples
#############################################################

ps <- prune_samples(
    sample_sums(ps) > 0,
    ps
)


#############################################################
## Check remaining samples
#############################################################

if (nsamples(ps) < 2) {
    
    stop(
        "Less than two samples remain after filtering."
    )
    
}

cat(
    "\nNumber of samples:",
    nsamples(ps),
    "\n"
)


#############################################################
## Convert phyloseq to SummarizedExperiment
#############################################################

se <- phyloseq_to_se(
    ps,
    tax_level = tax_level
)


#############################################################
## Convert counts to relative abundance
#############################################################

se.rel <- lefser::relativeAb(
    se
)


#############################################################
## Check relative abundance
#############################################################

abundance_sums <- colSums(
    SummarizedExperiment::assay(se.rel)
)

if (!all(
    abs(abundance_sums - 1) < 1e-8
)) {
    
    stop(
        "Relative abundance conversion failed: sample sums are not equal to 1."
    )
    
}


#############################################################
## Check groups
#############################################################

group_values <- as.character(
    SummarizedExperiment::colData(se.rel)[[group_var]]
)

group_values <- group_values[
    !is.na(group_values)
]

if (length(unique(group_values)) < 2) {
    
    stop(
        "LEfSe requires at least two groups."
    )
    
}


#############################################################
## Convert group to factor
#############################################################

SummarizedExperiment::colData(
    se.rel
)[[group_var]] <- factor(
    SummarizedExperiment::colData(
        se.rel
    )[[group_var]]
)


#############################################################
## Convert block to factor
#############################################################

if (!is.null(block_var)) {
    
    SummarizedExperiment::colData(
        se.rel
    )[[block_var]] <- factor(
        SummarizedExperiment::colData(
            se.rel
        )[[block_var]]
    )
    
}


#############################################################
## Display group sizes
#############################################################

cat(
    "\nGroup sizes:\n"
)

print(
    table(
        SummarizedExperiment::colData(
            se.rel
        )[[group_var]]
    )
)


#############################################################
## Display block sizes
#############################################################

if (!is.null(block_var)) {
    
    cat(
        "\nBlock sizes:\n"
    )
    
    print(
        table(
            SummarizedExperiment::colData(
                se.rel
            )[[block_var]]
        )
    )
    
}


#############################################################
## Run LEfSe
#############################################################

cat(
    "\nRunning LEfSe...\n"
)

lefse_result <- tryCatch(
    
    {
        
        if (is.null(block_var)) {
            
            lefser(
                se.rel,
                groupCol = group_var,
                lda.threshold = lda_cutoff,
                kruskal.threshold = kw_cutoff,
                checkAbundances = TRUE
            )
            
        } else {
            
            lefser(
                se.rel,
                groupCol = group_var,
                blockCol = block_var,
                lda.threshold = lda_cutoff,
                kruskal.threshold = kw_cutoff,
                wilcox.threshold = wilcox_cutoff,
                checkAbundances = TRUE
            )
            
        }
        
    },
    
    error = function(e) {
        
        message(
            "\nLEfSe failed:\n",
            e$message
        )
        
        return(NULL)
        
    }
    
)


#############################################################
## Check LEfSe results
#############################################################

if (is.null(lefse_result)) {
    
    warning(
        "LEfSe did not produce a result."
    )
    
} else {
    
    cat(
        "\nNumber of LEfSe biomarkers:",
        nrow(lefse_result),
        "\n"
    )
    
    if (nrow(lefse_result) > 0) {
        
        print(
            head(
                lefse_result,
                20
            )
        )
        
    } else {
        
        message(
            "\nNo taxa passed the selected LEfSe thresholds."
        )
        
    }
    
}


#############################################################
## Create output directory
#############################################################

if (!dir.exists(output_directory)) {
    
    dir.create(
        output_directory,
        recursive = TRUE
    )
    
}


#############################################################
## Export LEfSe results
#############################################################

if (!is.null(lefse_result)) {
    
    output_file <- file.path(
        output_directory,
        "LEfSe_results.csv"
    )
    
    write.csv(
        lefse_result,
        output_file,
        row.names = FALSE
    )
    
    cat(
        "\nResults saved to:\n",
        output_file,
        "\n"
    )
    
}


#############################################################
## Export analysis metadata
#############################################################

analysis_info <- data.frame(
    
    parameter = c(
        "tax_level",
        "group_variable",
        "block_variable",
        "LDA_threshold",
        "Kruskal_threshold",
        "Wilcoxon_threshold",
        "number_samples",
        "number_taxa"
    ),
    
    value = c(
        tax_level,
        group_var,
        ifelse(
            is.null(block_var),
            NA,
            block_var
        ),
        lda_cutoff,
        kw_cutoff,
        wilcox_cutoff,
        nsamples(se.rel),
        nrow(se.rel)
    )
    
)


write.csv(
    analysis_info,
    file.path(
        output_directory,
        "LEfSe_analysis_parameters.csv"
    ),
    row.names = FALSE
)


#############################################################
## Export group information
#############################################################

group_table <- as.data.frame(
    table(
        SummarizedExperiment::colData(
            se.rel
        )[[group_var]]
    )
)

colnames(group_table) <- c(
    "Group",
    "N"
)


write.csv(
    group_table,
    file.path(
        output_directory,
        "LEfSe_group_sizes.csv"
    ),
    row.names = FALSE
)


#############################################################
## End
#############################################################

cat(
    "\n=====================================================\n"
)

cat(
    "LEfSe analysis completed.\n"
)

cat(
    "=====================================================\n"
)
