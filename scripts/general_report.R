#! /usr/bin/env Rscript


##########################################
## LOAD LIBRARIES
##########################################
# Obtain this script directory
full.fpath <- normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                commandArgs())], '='))[2])

main_path_script <- dirname(full.fpath)
root_path <- file.path(main_path_script)
template_path <- file.path(root_path, "..", "templates")
# Load custom libraries
# devtools::load_all(file.path(root_path))

source_folder <- file.path(root_path, 'lib')
library(optparse)
library(Seurat)
source(file.path(source_folder, "preprocessing_library.R"))


##########################################
## OPTPARSE
##########################################

option_list <- list(
  optparse::make_option(c("-i", "--input"), type = "character",
              help="Input file with samples names"),
  optparse::make_option(c("-o", "--output"), type = "character",
              help="Folder where the report is written"),
  optparse::make_option(c("--filter"), type = "character",
              help="TRUE for using only detected cell-associated barcodes, FALSE for using all detected barcodes"),
  optparse::make_option(c("--mincells"), type = "integer",
              help="Min number of cells for which a feature was recorded"),
  optparse::make_option(c("--minfeats"), type = "integer",
              help="Min number of features for which a cell was recorded"),
  optparse::make_option(c("--minqcfeats"), type = "integer",
              help="Min number of features for which a cell was selected in QC"),
  optparse::make_option(c("--percentmt"), type = "integer",
              help="Max percentage of reads mapped to mitochondrial genes for which a cell is recorded"),
  optparse::make_option(c("--normalmethod"), type = "character",
              help="Method for normalization. LogNormalize, CLR or RC"),
  optparse::make_option(c("--scalefactor"), type = "integer",
              help="Scale factor for cell-level normalization"),
  optparse::make_option(c("--hvgs"), type = "integer",
              help="Number of HVG to be selected"),
  optparse::make_option(c("--ndims"), type = "integer",
              help="Number of PC to be used for clustering / UMAP / tSNE"),
  optparse::make_option(c("--dimheatmapcells"), type = "integer",
              help="Heatmap plots the 'extreme' cells on both ends of the spectrum"),
  optparse::make_option(c("--experiment_name"), type = "character",
              help="Experiment name"),
  optparse::make_option(c("--results_folder"), type = "character",
              help="Folder with the preprocessing results"),
  optparse::make_option(c("--resolution"), type = "double",
              help="Granularity of the clustering"),
  optparse::make_option(c("--integrative_analysis"), type = "character",
              help="TRUE if we want to integrate samples, FALSE otherwise"),
  optparse::make_option(c("--int_sec_cond"), type = "character",
              help="Secondary condition for splitted / grouped UMAPs, FALSE if not desired")
)  

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))


##########################################
## MAIN
##########################################
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
                           template = file.path(template_path,
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