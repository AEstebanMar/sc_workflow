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
options(future.globals.maxSize = 18000 * 1024^2)

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
  optparse::make_option(c("-p", "--project_name"), type = "character", default = NULL,
              help = "Project name"),
  optparse::make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output folder"),
  optparse::make_option("--filter", type = "character", default = NULL,
              help = "TRUE for using only detected cell-associated barcodes, FALSE for using all detected barcodes"),
  optparse::make_option("--mincells", type = "integer", default = NULL,
              help = "Min number of cells for which a feature was recorded"),
  optparse::make_option("--minfeats", type = "integer", default = NULL,
              help = "Min number of features for which a cell was recorded"),
  optparse::make_option("--minqcfeats", type = "integer", default = NULL,
              help = "Min number of features for which a cell was selected in QC"),
  optparse::make_option("--percentmt", type = "integer", default = NULL,
              help = "Max percentage of reads mapped to mitochondrial genes for which a cell is recorded"),
  optparse::make_option("--normalmethod", type = "character", default = NULL,
              help = "Method for normalization. LogNormalize, CLR or RC"),
  optparse::make_option("--scalefactor", type = "integer", default = NULL,
              help = "Scale factor for cell-level normalization"),
  optparse::make_option("--hvgs", type = "integer", default = NULL,
              help = "Number of HVG to be selected"),
  optparse::make_option("--ndims", type = "integer", default = NULL,
              help = "Number of PC to be used for clustering / UMAP / tSNE"),
  optparse::make_option("--dimheatmapcells", type = "integer", default = NULL,
              help = "Heatmap plots the 'extreme' cells on both ends of the spectrum"),
  optparse::make_option("--report_folder", type = "character", default = NULL,
              help = "Folder where the report is written"),
  optparse::make_option("--experiment_name", type = "character", default = NULL,
              help = "Experiment name"),
  optparse::make_option("--resolution", type = "double", default = NULL,
              help = "Granularity of the clustering"),
  optparse::make_option(c("-d", "--exp_design"), type = "character", default = NULL,
              help = "Input file with the experiment design"),
  optparse::make_option("--count_path", type = "character", default = NULL,
            help = "Count results folder"),
  optparse::make_option("--suffix", type = "character", help = "Suffix to specific file"),
  optparse::make_option("--samples_to_integrate", type = "character", default = "",
            help = "Comma-separated list of samples to integrate, will integrate all
                    samples if not specified."),
  optparse::make_option("--int_columns", type = "character", default = "",
            help = "Comma-separated list of conditions by which to perform integration
                    analysis. If empty, all conditions will be analysed."),
  optparse::make_option("--save_raw", type = "logical", default = FALSE, action = "store_true",
            help = "Save seurat object before analysis."),
  optparse::make_option("--clusters_annotation", type = "character", default = "",
            help = "Clusters annotation file."),
  optparse::make_option("--target_genes", type = "character", default = "",
            help = "Path to target genes table, or comma-separated list of target genes"),
  optparse::make_option("--cpus", type = "double", default = 1,
            help = "Provided CPUs"),
  optparse::make_option("--imported_counts", type = "character", default = "",
            help = "Imported counts directory"),
  optparse::make_option("--DEG_columns", type = "character", default = "",
            help = "Columns for DEG analysis"),
  optparse::make_option("--cell_types_annotation", type = "character", default = "",
            help = "Columns for DEG analysis")
)  

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

plan("multicore", workers = opt$cpu)

##########################################
## MAIN
##########################################
if(opt$target_genes == ""){
  warning("No target genes provided")
  target_genes <- NULL
} else if(file.exists(opt$target_genes)) {
  target_genes <- read_and_format_markers(opt$target_genes)
} else {
  target_genes <- list(Custom = strsplit(opt$target_genes, split = ",")[[1]])
}

exp_design <- read.table(opt$exp_design, sep = "\t", header = TRUE)
if(is.na(opt$int_columns)) {
  warning("No conditions specified for integration. Analysing every condition")
  int_columns <- colnames(exp_design)[!colnames(exp_design)=="sample"]
} else {
  int_columns <- unlist(strsplit(opt$int_columns, ","))
}

message(paste0("Selected ", length(int_columns), " condition(s) for analysis: ", paste0(int_columns, collapse = ", ")))

# Input reading and integration variables setup
if(opt$samples_to_integrate == "") {
  samples <- exp_design$sample
} else {
  samples <- strsplit(opt$samples_to_integrate, ",")[[1]]
}

dimreds_to_do <- c("pca") # For dimensionality reduction
embeds <- "harmony"

out_path = file.path(opt$report_folder, opt$experiment_name)

# if(opt$imported_counts == "") {
#   merged_seu <- merge_seurat(project = opt$project_name, samples = samples, exp_design = exp_design,
#                             suffix = opt$suffix, count_path = opt$count_path)  
# } else {
#   merged_seu <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(opt$imported_counts, gene.column = 1),
#                                            project = opt$experiment_name, min.cells = 1, min.features = 1)
#   merged_seu_meta <- read.table(file.path(opt$imported_counts, "meta.tsv"), sep = "\t", header = TRUE)
#   rownames(merged_seu_meta) <- colnames(merged_seu)
#   merged_seu <- AddMetaData(merged_seu, merged_seu_meta, row.names("Cell_ID"))
# }

# message('Normalizing data')
# merged_seu <- Seurat::NormalizeData(merged_seu, normalization.method = opt$normalmethod,
#                        scale.factor = opt$scalefactor, verbose = FALSE)
# message('Scaling data')
# merged_seu <- Seurat::ScaleData(merged_seu, features = rownames(merged_seu), verbose = FALSE)

# if(opt$save_raw) {
#     saveRDS(seu, file = file.path(out_path, paste0(opt$experiment_name, ".before.seu.RDS")))
#   }

# message('Analyzing full experiment')

# global_seu <- analyze_seurat(raw_seu = merged_seu, out_path = out_path, minqcfeats = opt$minqcfeats,
#                              percentmt = opt$percentmt, hvgs = opt$hvgs, ndims = opt$ndims,
#                              resolution = opt$resolution, dimreds_to_do = dimreds_to_do,
#                              embeds = embeds, integrate = TRUE)

# message("--------------------------------------------")
# message("---------Writing global report---------")
# message("--------------------------------------------")

# write_seurat_report(all_seu = global_seu, percentmt = opt$percentmt, template = file.path(template_path,
#                     "preprocessing_report.Rmd"), out_path = out_path, minqcfeats = opt$minqcfeats,
#                     intermediate_files = "int_files", hvgs = opt$hvgs, resolution = opt$resolution)

# message('Starting integration analysis')

# message("Integrating seurat object")

# int_seu <- integrate_seurat(seu = merged_seu, clusters_annotation = opt$clusters_annotation,
#                             ndims = opt$ndims, resolution = opt$resolution, embeds = embeds,
#                             minqcfeats = opt$minqcfeats, percentmt = opt$percentmt, hvgs = opt$hvgs,
#                             scalefactor = opt$scalefactor, normalmethod = opt$normalmethod)

# saveRDS(int_seu, "int_seu.rds")
# stop('test')
int_seu <- readRDS('mouse_int_seu.rds')

DEG_list <- list()

if(opt$DEG_columns == "") {
  DEG_conditions <- int_columns
} else {
  DEG_conditions <- unlist(strsplit(opt$DEG_columns, split = ","))
}
for(condition in DEG_conditions) {
  message(paste0("Calculating DEGs for condition ", condition))
  condition_DEGs <- get_sc_markers(seu = int_seu, cond = condition, DEG = TRUE)
  DEG_list[[condition]] <- condition_DEGs
}

saveRDS(DEG_list, "mouse_DEG_list.rds")

conserved_markers <- list()
markers_list <- list()
if(length(int_columns) == 1) {
  markers <- get_sc_markers(seu = int_seu, cond = int_columns, DEG = FALSE, top = 200)
  markers <- collapse_markers(markers)
} else {
  markers <- Seurat::FindAllMarkers(int_seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
}

saveRDS(markers, "mouse_markers.rds")

if(opt$annotated_cell_types != "" & opt$annotated_cell_types != "") {
  anno_table <- read.table(opt$annotated_cell_types, sep = "\t", header = TRUE)
  annotated_clusters <- match_cell_types(markers, anno_table)
  int_seu$seu <- annotate_clusters(int_seu$seu, annotated_clusters$cell_types)
} else {
  message("No cell type data provided. Clusters cannot be annotated")
}

message("--------------------------------------------")
message("---------Writing integration report---------")
message("--------------------------------------------")
write_integration_report(int_seu = int_seu, template_folder = template_path, DEG_list = DEG_list,
                         markers_list = markers_list, output_dir = opt$report_folder, source_folder = source_folder,
                         target_genes = target_genes, name = opt$project_name, int_columns = int_columns)
message(paste0("Report written in ", opt$report_folder))
