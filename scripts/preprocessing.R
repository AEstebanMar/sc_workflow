#! /usr/bin/env Rscript


# Sergio Alías, 20230606
# Last modified 20231226


#################################
###   STAGE 3 PREPROCESSING   ###
###   preprocessing.R         ###
#################################

#################
### Libraries ###
#################

library(optparse)
library(Seurat)
library(scCustomize)

###################
### Custom libs ###
###################

root_path <- Sys.getenv("CODE_PATH") # daemon (TODO decide what to to with this)
source(file.path(root_path, "R", "preprocessing_library.R"))

############
### Args ###
############

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help="Input folder with 10X data"),
  make_option(c("-o", "--output"), type = "character",
              help="Output folder"),
  make_option(c("-n", "--name"), type = "character",
              help="Sample name"),
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
  make_option(c("--report_folder"), type = "character",
              help="Folder where the report is written"),
  make_option(c("--experiment_name"), type = "character",
              help="Experiment name"),
  make_option(c("--resolution"), type = "double",
              help="Granularity of the clustering"),
  make_option(c("--integrate"), type = "logical", default=FALSE, action = "store_true",
              help="Perform integrative analysis")
)  


opt <- parse_args(OptionParser(option_list = option_list))


############
### Main ###
############

out_path = file.path(opt$report_folder, paste0(opt$experiment_name, '.', opt$name))
# Input selection

# Input reading and integration variables setup

if (integrate) {
  seu <- readRDS(input)
  dimreds_to_do <- c("pca") # For dimensionality reduction
  embeddings_to_use <- "harmony"
} else {
  input <- file.path(opt$input, ifelse(opt$filter, 
      "filtered_feature_bc_matrix", "raw_feature_bc_matrix"))
  seu <- read_input(name = opt$name, 
                    input = input,
                    mincells = opt$mincells,
                    minfeats = opt$minfeats)
  dimreds_to_do <- c("pca", "tsne", "umap") # For dimensionality reduction
  embeddings_to_use <- "pca"
}

all_seu <- main_preprocessing_analysis(seu = seu,
                            out_path = out_path,
                            minqcfeats = opt$minqcfeats,
                            percentmt = opt$percentmt,
                            normalmethod = opt$normalmethod,
                            scalefactor = opt$scalefactor,
                            hvgs = opt$hvgs,
                            ndims = opt$ndims,
                            resolution = opt$resolution,
                            dimreds_to_do = dimreds_to_do,
                            embeddings_to_use = embeddings_to_use,
                            integrate = as.logical(opt$integrate))

write_preprocessing_report(all_seu = all_seu,
                           template = file.path(root_path, 
                                                "templates",
                                                "preprocessing_report.Rmd"),
                           out_path = out_path,
                           intermediate_files = "int_files",
                           minqcfeats = opt$minqcfeats,
                           percentmt = opt$percentmt,
                           hvgs = opt$hvgs,
                           resolution = opt$resolution)