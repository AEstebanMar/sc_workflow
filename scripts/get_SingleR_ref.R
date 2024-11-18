#! /usr/bin/env Rscript


option_list <- list(
  optparse::make_option(c("-r", "--reference"), type = "character",
            help = "Celldex reference to use."),
  optparse::make_option(c("-v", "--version"), type = "character",
            help = "Celldex version of reference."),
  optparse::make_option(c("-o", "--output"), type = "character",
            help = "Output directory."),
  optparse::make_option("--replace", type = "logical", default = FALSE, action = "store_true",
  			help = "Replace existing directory"),
  optparse::make_option("--verbose", type = "logical", default = TRUE, action = "store_true",
  	help = "Display progress"),
  optparse::make_option("--quiet", type = "logical", default = FALSE, action = "store_false",
  	dest = "verbose", help = "Display progress"),
  optparse::make_option("--database", type = "character", help = "Database to consult
    (\"celldex\" or \"scRNAseq\") or a path to local reference"),
  optparse::make_option("--ref_label", type = "character", default = NULL,
    help = "Column of reference metadata to use for annotation. Only used in scRNAseq mode"),
  optparse::make_option("--cpu", type = "integer", default = NULL,
    help = "CPUs to use when calling data.table::fread")
)  

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
output <- file.path(opt$output, paste(opt$reference))
if(opt$version != "") {
  output <- paste(output, opt$version, sep = "_")
}

if(is.null(opt$reference)) {
	stop('No reference provided. Please see get_celldex_ref.R --help')
}
if(is.null(opt$output)) {
	stop('No output path provided. Please see get_celldex_ref.R --help')
}
if(file.exists(opt$database)) {
  metadata <- read.table(file.path(opt$database, "metadata/meta*.txt"), sep = '\t')
  expr_files <- Sys.glob(paste0(opt$database, "/expression/*txt*"))
  if(length(expr_files) > 0) {
    data.table::setDTthreads(threads = opt$CPU)
    expr_matrices <- vector(mode = "list", length = length(expr_files))
    for(file in seq(expr_files)) {
      message(paste0("Loading expression file ", file, " of ", length(expr_files)))
      expr_matrix <- as.matrix(data.table::fread(expr_files[file], verbose = TRUE))
      message("File loaded. Processing and converting to sparse matrix.")
      row_names <- expr_matrix[, 1]
      expr_matrix <- expr_matrix[, -1]
      expr_matrix <- apply(expr_matrix, 2, as.numeric)
      rownames(expr_matrix) <- row_names
      sparse_matrix <- Matrix::Matrix(expr_matrix, sparse = TRUE)
      sparse_matrix <- sparse_matrix[, Matrix::colSums(sparse_matrix) > 0]
      expr_matrices[[file]] <- sparse_matrix
      message(paste0("File ", file, " converted"))
    }
    read_sparse_matrix <- do.call(cbind, expr_matrices)
  } else {
    read_sparse_matrix <- NULL
  }
  tenX_dirs <- dirname(Sys.glob("./expression/*/*mtx*"))
  if(length(tenX_dirs) > 0) {
    tenX_matrices <- vector(mode = "list", length = length(tenX_dirs))
    for(tenX_dir in seq(tenX_dirs)) {
      message(paste0("Loading 10X directory ", tenX_dir, " of ", length(tenX_dirs)))
      tenX_matrix <- Seurat::Read10X(tenX_dirs[tenX_dir], gene.column = 1)
      tenX_matrices[[tenX_dir]] <- tenX_matrix
    }
    message("10X matrices loaded. Merging (this may take a while)")
    merged_tenX_matrix <- SeuratObject::RowMergeSparseMatrices(tenX_matrices[[1]], tenX_matrices[-1])
  } else {
    read_sparse_matrix <- NULL
  }
  save_image('ref_testing.RData')
  if(!is.null(read_sparse_matrix) & !is.null(merged_tenX_matrix)) {
   final_sparse_matrix <- SeuratObject::RowMergeSparseMatrices(read_sparse_matrix, merged_tenX_matrix)  
  } else {
    if(is.null(read_sparse_matrix)) {
      final_sparse_matrix <- merged_tenX_matrix
    } else {
      final_sparse_matrix <- read_sparse_matrix
    }
  }
  final_sparse_matrix <- final_sparse_matrix[, order(colnames(final_sparse_matrix))]
  metadata <- metadata[metadata$NAME %in% colnames(final_sparse_matrix), ]
  metadata <- metadata[order(metadata$NAME), ]
  rownames(metadata) <- metadata$NAME
  metadata <- metadata[, -1]
  ref <- SummarizedExperiment::SummarizedExperiment(assays=list(counts = final_sparse_matrix), colData = metadata)
  ref <- ref[, !is.na(ref[[opt$ref_label]]) & ref[[opt$ref_label]] != "unclear"]
  ref <- scater::logNormCounts(ref)
}
if(opt$database == "scRNAseq") {
  ref <- scRNAseq::fetchDataset(opt$reference, opt$version)
  # Removing unlabelled cells or cells without a clear label.
  ref <- ref[, !is.na(ref[[opt$ref_label]]) & ref[[opt$ref_label]]!="unclear"] 
  ref <- scater::logNormCounts(ref)
}
if(opt$database == "celldex") {
  ref <- celldex::fetchReference(opt$reference, opt$version)
}

HDF5Array::saveHDF5SummarizedExperiment(x = ref, dir = output, verbose = opt$verbose,
										replace = opt$replace)

message(paste0("Reference saved successfully in ", output))
