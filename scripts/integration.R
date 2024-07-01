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

library(future)
options(future.globals.maxSize = 6000 * 1024^2)

sc_source_folder <- file.path(root_path, 'lib')
library(Seurat)
library(scCustomize)
source(file.path(sc_source_folder, "preprocessing_library.R"))

# Temporary path until we have htmlreportR installed
devtools::load_all('~/dev_R/htmlreportR')
source_folder <- file.path(find.package("htmlreportR"), "inst")

##########################################
## OPTPARSE
##########################################

option_list <- list(
  optparse::make_option(c("-p", "--project_name"), type = "character",
              help = "Project name"),
  optparse::make_option(c("-o", "--output"), type = "character",
              help = "Output folder"),
  optparse::make_option("--filter", type = "character",
              help = "TRUE for using only detected cell-associated barcodes, FALSE for using all detected barcodes"),
  optparse::make_option("--mincells", type = "integer",
              help = "Min number of cells for which a feature was recorded"),
  optparse::make_option("--minfeats", type = "integer",
              help = "Min number of features for which a cell was recorded"),
  optparse::make_option("--minqcfeats", type = "integer",
              help = "Min number of features for which a cell was selected in QC"),
  optparse::make_option("--percentmt", type = "integer",
              help = "Max percentage of reads mapped to mitochondrial genes for which a cell is recorded"),
  optparse::make_option("--normalmethod", type = "character",
              help = "Method for normalization. LogNormalize, CLR or RC"),
  optparse::make_option("--scalefactor", type = "integer",
              help = "Scale factor for cell-level normalization"),
  optparse::make_option("--hvgs", type = "integer",
              help = "Number of HVG to be selected"),
  optparse::make_option("--ndims", type = "integer",
              help = "Number of PC to be used for clustering / UMAP / tSNE"),
  optparse::make_option("--dimheatmapcells", type = "integer",
              help = "Heatmap plots the 'extreme' cells on both ends of the spectrum"),
  optparse::make_option("--report_folder", type = "character",
              help = "Folder where the report is written"),
  optparse::make_option("--experiment_name", type = "character",
              help = "Experiment name"),
  optparse::make_option("--resolution", type = "double",
              help = "Granularity of the clustering"),
  optparse::make_option(c("-d", "--exp_design"), type = "character",
              help = "Input file with the experiment design"),
  optparse::make_option("--count_path", type = "character",
            help = "Count results folder"),
  optparse::make_option("--suffix", type = "character", help = "Suffix to specific file"),
  optparse::make_option("--int_columns", type = "character", default = NULL,
            help = "Columns by which seurat object will be subsetted in integration
                    analysis"),
  optparse::make_option("--samples_to_integrate", type = "character", default = NULL,
            help = "Comma-separated list of samples to integrate, will integrate all
                    samples if not specified."),
  optparse::make_option("--save_raw", type = "logical", default = FALSE, action = "store_true",
            help = "Save seurat object before analysis."),
  optparse::make_option("--annotation_dir", type = "character", default = NULL,
            help = "Path to directory containing cluster annotation files."),
  optparse::make_option("--target_genes", type = "character", default = NULL,
            help = "Path to target genes table."),
  optparse::make_option("--cpus", type = "double", default = 1,
            help = "Provided CPUs")
)  

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

plan("multicore", workers = opt$cpu)

##########################################
## MAIN
##########################################

target_genes <- read_and_format_markers(opt$target_genes)
exp_design <- read.table(opt$exp_design, sep = "\t", header = TRUE)
subset_columns <- strsplit(opt$int_columns, ",")[[1]]

# Input reading and integration variables setup
if(opt$samples_to_integrate == "") {
  samples <- exp_design$code
} else {
  samples <- strsplit(opt$samples_to_integrate, ",")[[1]]
}

dimreds_to_do <- c("pca") # For dimensionality reduction
embeddings_to_use <- "harmony"

out_path = file.path(opt$report_folder, opt$experiment_name)

merged_seu <- merge_seurat(project = opt$project_name, samples = samples, exp_design = exp_design,
                            suffix = opt$suffix, count_path = opt$count_path)
message('Normalizing data')
merged_seu <- NormalizeData(merged_seu, normalization.method = opt$normalmethod,
                       scale.factor = opt$scalefactor, verbose = FALSE)
message('Scaling data')
merged_seu <- ScaleData(merged_seu, features = rownames(merged_seu), verbose = FALSE)

if(opt$save_raw) {
    saveRDS(seu, file = file.path(out_path, paste0(opt$experiment_name, ".before.seu.RDS")))
  }

message('Analyzing full experiment')

global_seu <- analyze_seurat(raw_seu = merged_seu, out_path = out_path, minqcfeats = opt$minqcfeats,
                             percentmt = opt$percentmt, hvgs = opt$hvgs, ndims = opt$ndims,
                             resolution = opt$resolution, dimreds_to_do = dimreds_to_do,
                             embeddings_to_use = embeddings_to_use, integrate = TRUE)

message("--------------------------------------------")
message("---------Writing global report---------")
message("--------------------------------------------")

write_seurat_report(all_seu = global_seu, percentmt = opt$percentmt, template = file.path(template_path,
                    "preprocessing_report.Rmd"), out_path = out_path, minqcfeats = opt$minqcfeats,
                    intermediate_files = "int_files", hvgs = opt$hvgs, resolution = opt$resolution)

message('Starting integration analysis')

for(column in subset_columns) {
  message(paste0("Integrating variable \"", column, "\""))
  sec_column <- subset_columns[subset_columns != column]
  if(length(sec_column) == 0) {
    sec_column <- NULL
  }
  # comparison <- readRDS('comparison.rds')
  comparison <- compare_subsets(seu = merged_seu, annotation_dir = opt$annotation_dir, subset_column = column,
                                exp_design = exp_design, ndims = opt$ndims, resolution = opt$resolution,
                                embeddings_to_use = embeddings_to_use, minqcfeats = opt$minqcfeats,
                                percentmt = opt$percentmt, hvgs = opt$hvgs, scalefactor = opt$scalefactor,
                                normalmethod = opt$normalmethod)
  message("--------------------------------------------")
  message("---------Writing integration report---------")
  message("--------------------------------------------")
  write_integration_report(comparison = comparison, template_folder = template_path,
                           output_dir = opt$report_folder, source_folder = source_folder,
                           target_genes = target_genes, name = column, int_column = column,
                           sec_column = sec_column)
  message(paste0("Report written in ", opt$report_folder))
}
