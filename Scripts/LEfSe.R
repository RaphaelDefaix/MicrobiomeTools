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



LEfSe_results.csv
