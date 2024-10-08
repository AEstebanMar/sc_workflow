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
options(future.globals.maxSize = 18 * 1024^6)

sc_source_folder <- file.path(root_path, 'lib')
source(file.path(sc_source_folder, "sc_library.R"))
source(file.path(sc_source_folder, "main_analyze_seurat.R"))

# Temporary path until we have htmlreportR installed
devtools::load_all('~/dev_R/htmlreportR')
source_folder <- file.path(find.package("htmlreportR"), "inst")

##########################################
## OPTPARSE
##########################################

option_list <- list(
  optparse::make_option(c("-p", "--project_name"), type = "character", default = NULL,
              help = "Project name."),
  optparse::make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output folder."),
  optparse::make_option("--filter", type = "character", default = NULL,
              help = "TRUE for using only detected cell-associated barcodes,
                      FALSE for using all detected barcodes."),
  optparse::make_option("--mincells", type = "integer", default = NULL,
              help = "Min number of cells for which a feature was recorded."),
  optparse::make_option("--minfeats", type = "integer", default = NULL,
              help = "Min number of features for which a cell was recorded."),
  optparse::make_option("--minqcfeats", type = "integer", default = NULL,
              help = "Min number of features for which a cell was selected in QC."),
  optparse::make_option("--percentmt", type = "integer", default = NULL,
              help = "Max percentage of reads mapped to mitochondrial genes for which a cell is recorded."),
  optparse::make_option("--normalmethod", type = "character", default = NULL,
              help = "Method for normalization. LogNormalize, CLR or RC."),
  optparse::make_option("--scalefactor", type = "integer", default = NULL,
              help = "Scale factor for cell-level normalization."),
  optparse::make_option("--hvgs", type = "integer", default = NULL,
              help = "Number of HVG to be selected."),
  optparse::make_option("--ndims", type = "integer", default = NULL,
              help = "Number of PC to be used for clustering / UMAP / tSNE"),
  optparse::make_option("--dimheatmapcells", type = "integer", default = NULL,
              help = "Heatmap plots the 'extreme' cells on both ends of the spectrum."),
  optparse::make_option("--resolution", type = "double", default = NULL,
              help = "Granularity of the clustering."),
  optparse::make_option(c("-d", "--exp_design"), type = "character", default = NULL,
              help = "Input file with the experiment design."),
  optparse::make_option("--count_path", type = "character", default = NULL,
            help = "Count results folder."),
  optparse::make_option("--suffix", type = "character", help = "Suffix to specific file"),
  optparse::make_option("--samples_to_integrate", type = "character", default = "",
            help = "Path to file containing samples to be processed."),
  optparse::make_option("--int_columns", type = "character", default = "",
            help = "Comma-separated list of conditions by which to perform integration
                    analysis. If empty, all conditions will be analysed."),
  optparse::make_option("--cluster_annotation", type = "character", default = "",
            help = "Clusters annotation file."),
  optparse::make_option("--target_genes", type = "character", default = "",
            help = "Path to target genes table, or comma-separated list of target genes."),
  optparse::make_option("--cpu", type = "double", default = 1,
            help = "Provided CPUs."),
  optparse::make_option("--imported_counts", type = "character", default = "",
            help = "Imported counts directory."),
  optparse::make_option("--DEG_columns", type = "character", default = "",
            help = "Columns for DEG analysis."),
  optparse::make_option("--cell_annotation", type = "character", default = "",
            help = "Cell types annotation file. Will be used to dynamically
                    annotate clusters."),
  optparse::make_option("--p_adj_cutoff", type = "numeric", default = "5e-3",
            help = "Adjusted p-value cutoff."),
  optparse::make_option("--verbose", type = "logical", default = FALSE, action = "store_true",
            help = "Verbosity of base Seurat and harmony function calls."),
  optparse::make_option("--reduce", type = "logical", default = FALSE, action = "store_true",
            help = "Randomly subset seurat object to 3000 cells, for quick testing.")
)  

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

plan("multicore", workers = opt$cpu)

##########################################
## MAIN
##########################################
if(opt$cluster_annotation != "") {
  cluster_annotation <- read.table(opt$cluster_annotation, sep = "\t", header = TRUE)
} else {
  cluster_annotation <- NULL
}
if(opt$cell_annotation != "") {
  cell_annotation <- read.table(opt$cell_annotation, sep = "\t", header = TRUE)
} else {
  cell_annotation <- NULL
}

if(opt$target_genes == ""){
  warning("No target genes provided")
  target_genes <- NULL
} else if(file.exists(opt$target_genes)) {
  target_genes <- read_and_format_targets(opt$target_genes)
} else {
  target_genes <- list(Custom = strsplit(opt$target_genes, split = ";")[[1]])
}

exp_design <- read.table(opt$exp_design, sep = "\t", header = TRUE)

if(opt$int_columns == "") {
  warning("No conditions specified for integration. Analysing every condition")
  int_columns <- colnames(exp_design)[!colnames(exp_design)=="sample"]
} else {
  int_columns <- unlist(strsplit(opt$int_columns, ","))
}

if(opt$DEG_columns == "") {
  DEG_columns <- int_columns
} else {
  DEG_columns <- opt$DEG_columns
}

message(paste0("Selected ", length(int_columns), " condition(s) for analysis: ", paste0(int_columns, collapse = ", ")))

# Input reading and integration variables setup
if(opt$samples_to_integrate == "") {
  samples <- exp_design$sample
} else {
  samples <- read.table(opt$samples_to_integrate, sep = "\t", header = FALSE)[[1]]
}

exp_design <- exp_design[exp_design$sample %in% samples, ]

if(opt$imported_counts == "") {
  merged_seu <- merge_seurat(project = opt$project_name, exp_design = exp_design,
                             suffix = opt$suffix, count_path = opt$count_path)  
} else {
  merged_seu <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(opt$imported_counts, gene.column = 1),
                                           project = opt$project_name, min.cells = 1, min.features = 1)
  merged_seu_meta <- read.table(file.path(opt$imported_counts, "meta.tsv"), sep = "\t", header = TRUE)
  rownames(merged_seu_meta) <- colnames(merged_seu)
  merged_seu <- Seurat::AddMetaData(merged_seu, merged_seu_meta, row.names("Cell_ID"))
}

if(opt$reduce) {
  message('Downsampling seurat object')
  merged_seu <- downsample_seurat(merged_seu, cells = 500, features = 5000, keep = unlist(target_genes))
}

message("Analyzing seurat object")

final_results <- main_analyze_seurat(seu = merged_seu, cluster_annotation = cluster_annotation,
                                     ndims = opt$ndims, resolution = opt$resolution, int_columns = int_columns,
                                     cell_annotation = cell_annotation, DEG_columns = DEG_columns,
                                     minqcfeats = opt$minqcfeats, percentmt = opt$percentmt, hvgs = opt$hvgs,
                                     scalefactor = opt$scalefactor, normalmethod = opt$normalmethod,
                                     p_adj_cutoff = opt$p_adj_cutoff, verbose = opt$verbose, sigfig = 2,
                                     output = opt$output, integrate = TRUE, query = unlist(target_genes),
                                     reduce = opt$reduce, save_RDS = TRUE)

message("-----------------------------------")
message("---------Writing QC report---------")
message("-----------------------------------")

write_seurat_report(final_results = final_results, template_folder = template_path,
                    template = "qc_template.txt", output = file.path(opt$output, "report"),
                    source_folder = source_folder, target_genes = target_genes,
                    name = opt$project_name, out_name = "qc_report.html", use_canvas = TRUE)

message("--------------------------------------------")
message("---------Writing integrative report---------")
message("--------------------------------------------")

write_seurat_report(final_results = final_results, template_folder = template_path,
                    output = file.path(opt$output, "report"), source_folder = source_folder,
                    target_genes = target_genes, name = opt$project_name,
                    int_columns = int_columns, cell_annotation = cell_annotation,
                    template = "analysis_template.txt", out_name = "integration_report.html")
