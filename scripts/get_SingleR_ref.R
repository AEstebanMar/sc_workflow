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
  optparse::make_option("--mode", type = "character", help = "Database to consult
    (\"celldex\" or \"scRNAseq\") or create from local files (\"custom\")"),
  optparse::make_option("--ref_label", type = "character", default = NULL,
    help = "Column of reference metadata to use for annotation. Only used in scRNAseq mode")
)  

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

output <- file.path(opt$output, paste(opt$reference, opt$version, sep = "_"))

if(is.null(opt$reference)) {
	stop('No reference provided. Please see get_celldex_ref.R --help')
}
if(is.null(opt$output)) {
	stop('No output path provided. Please see get_celldex_ref.R --help')
}
if(opt$mode == "custom") {
  # Not yet functional. Testing in create_SCP_ref.R. Once it is working, it will
  # be moved to this block
  custom_things
}
if(opt$mode == "scRNAseq") {
  ref <- scRNAseq::fetchDataset(opt$reference, opt$version)
  # Removing unlabelled cells or cells without a clear label.
  ref <- ref[, !is.na(ref[[opt$ref_label]]) & ref[[opt$ref_label]]!="unclear"] 
  ref <- scater::logNormCounts(ref)
}
if(opt$mode == "celldex") {
  ref <- celldex::fetchReference(opt$reference, opt$version)
}

HDF5Array::saveHDF5SummarizedExperiment(x = ref, dir = output, verbose = opt$verbose,
										replace = opt$replace)
