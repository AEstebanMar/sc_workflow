#! /usr/bin/env Rscript

# Sergio Alías, 20230627
# Last modified 20231216


################################################
###   STAGE 3 GENERAL PREPROCESSING REPORT   ###
###   general_report.R                       ###
################################################

#################
### Libraries ###
#################

library(optparse)
library(Seurat)


###################
### Custom libs ###
###################

root_path <- Sys.getenv("CODE_PATH") # daemon
source(file.path(root_path, "R", "preprocessing_library.R"))

############
### Args ###
############

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help="Input file with samples names"),
  make_option(c("-o", "--output"), type = "character",
              help="Folder where the report is written"),
  make_option(c("--filter"), type = "character",
              help="TRUE for using only detected cell-associated barcodes, FALSE for using all detected barcodes"),
  make_option(c("--mincells"), type = "integer",
              help="Min number of cells for which a feature was recorded"),
  make_option(c("--minfeats"), type = "integer",
              help="Min number of features for which a cell was recorded"),
  make_option(c("--minqcfeats"), type = "integer",
              help="Min number of features for which a cell was selected in QC"),
  make_option(c("--percentmt"), type = "integer",
              help="Max percentage of reads mapped to mitochondrial genes for which a cell is recorded"),
  make_option(c("--normalmethod"), type = "character",
              help="Method for normalization. LogNormalize, CLR or RC"),
  make_option(c("--scalefactor"), type = "integer",
              help="Scale factor for cell-level normalization"),
  make_option(c("--hvgs"), type = "integer",
              help="Number of HVG to be selected"),
  make_option(c("--ndims"), type = "integer",
              help="Number of PC to be used for clustering / UMAP / tSNE"),
  make_option(c("--dimheatmapcells"), type = "integer",
              help="Heatmap plots the 'extreme' cells on both ends of the spectrum"),
  make_option(c("--experiment_name"), type = "character",
              help="Experiment name"),
  make_option(c("--results_folder"), type = "character",
              help="Folder with the preprocessing results"),
  make_option(c("--resolution"), type = "double",
              help="Granularity of the clustering"),
  make_option(c("--integrative_analysis"), type = "character",
              help="TRUE if we want to integrate samples, FALSE otherwise"),
  make_option(c("--int_sec_cond"), type = "character",
              help="Secondary condition for splitted / grouped UMAPs, FALSE if not desired")
)  

opt <- parse_args(OptionParser(option_list = option_list))


############
### Main ###
############
both_seu_paths <- Sys.glob(paste0(opt$results_folder, "seu.RDS"))
seu_paths <- grep("before", both_seu_paths, invert = TRUE, value = TRUE)
raw_seu_paths <- grep("before", both_seu_paths, value = TRUE)
marker_gene_paths <- Sys.glob(paste0(opt$results_folder, ".markers.RDS"))

seu <- lapply(seu_paths, readRDS)
before.seu <- lapply(raw_seu_paths, readRDS)
marker_gene_list <- lapply(marker_gene_paths, readRDS)

main_folder <- opt$results_folder #Reuse
experiment_name <- opt$experiment_name

if (isTRUE(opt$integrative_analysis)){
  report_name <- "All integrated samples"
} else {
  report_name <- "All samples"
}

write_preprocessing_report(name = report_name,
                           experiment = experiment_name,
                           template = file.path(root_path, 
                                                "templates",
                                                "preprocessing_report.Rmd"),
                           outdir = opt$output,
                           intermediate_files = "int_files",
                           minqcfeats = opt$minqcfeats,
                           percentmt = opt$percentmt,
                           hvgs = opt$hvgs,
                           resolution = opt$resolution,
                           all_seu = list(seu, before.seu, marker_gene_list)
                           #integrate = opt$integrative_analysis,
                           #sec_cond = opt$int_sec_cond
                           )