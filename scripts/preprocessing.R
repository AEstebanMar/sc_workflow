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
library(scCustomize)
source(file.path(source_folder, "main_analyze_seurat.R"))
source(file.path(source_folder, "preprocessing_library.R"))

##########################################
## OPTPARSE
##########################################

option_list <- list(
  optparse::make_option(c("-i", "--input"), type = "character",
              help="Input folder with 10X data"),
  optparse::make_option(c("-o", "--output"), type = "character",
              help="Output folder"),
  optparse::make_option(c("-n", "--name"), type = "character",
              help="Sample name"),
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
  optparse::make_option(c("--report_folder"), type = "character",
              help="Folder where the report is written"),
  optparse::make_option(c("--experiment_name"), type = "character",
              help="Experiment name"),
  optparse::make_option(c("--resolution"), type = "double",
              help="Granularity of the clustering"),
  optparse::make_option(c("--exp_design"), type = "character",
              help="Path to experiment design file"),
  optparse::make_option("--target_genes", type = "character", default = "",
              help = "Path to target genes table, or comma-separated list of target genes")
)  


opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

##########################################
## MAIN
##########################################

out_path = file.path(opt$report_folder, paste0(opt$experiment_name, '.', opt$name))
# Input selection

# Input reading and integration variables setup

input <- file.path(opt$input, ifelse(opt$filter, "filtered_feature_bc_matrix",
                                                 "raw_feature_bc_matrix"))

exp_design <- read.table(opt$exp_design, sep = "\t", header = TRUE)

if(opt$target_genes == ""){
  warning("No target genes provided")
  target_genes <- NULL
} else if(file.exists(opt$target_genes)) {
  target_genes <- read_and_format_targets(opt$target_genes)
} else {
  target_genes <- list(Custom = strsplit(opt$target_genes, split = ";")[[1]])
}

seu <- read_input(name = opt$name, input = input, mincells = opt$mincells,
                  minfeats = opt$minfeats, exp_design = exp_design)

all_seu <- main_analyze_seurat(seu = seu, output = opt$report_folder, hvgs = opt$hvgs,
                               minqcfeats = opt$minqcfeats, ndims = opt$ndims,
                               percentmt = opt$percentmt,
                               resolution = opt$resolution,
                               integrate = FALSE, query = unlist(target_genes))

write_sergio_report(all_seu = all_seu, minqcfeats = opt$minqcfeats,
                    template = file.path(template_path,
                                         "preprocessing_report.Rmd"),
                    out_path = out_path, resolution = opt$resolution,
                    intermediate_files = "int_files",
                    percentmt = opt$percentmt, hvgs = opt$hvgs)
