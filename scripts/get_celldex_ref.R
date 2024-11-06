#! /usr/bin/env Rscript


option_list <- list(
  optparse::make_option(c("-r", "--reference"), type = "character",
            help = "Celldex reference to use."),
  optparse::make_option(c("-v", "--version"), type = "character", default = "2024-02-26",
            help = "Celldex version of reference."),
  optparse::make_option(c("-o", "--output"), type = "character",
            help = "Output directory."),
  optparse::make_option("--replace", type = "logical", default = FALSE, action = "store_true",
  			help = "Replace existing directory"),
  optparse::make_option("--verbose", type = "logical", default = TRUE, action = "store_true",
  	help = "Display progress"),
  optparse::make_option("--quiet", type = "logical", default = FALSE, action = "store_false",
  	dest = "verbose", help = "Display progress")
)  

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if(is.null(opt$reference)) {
	stop('No reference provided. Please see get_celldex_ref.R --help')
}
if(is.null(opt$output)) {
	stop('No output path provided. Please see get_celldex_ref.R --help')
}

output <- file.path(opt$output, opt$reference)
ref <- celldex::fetchReference(opt$reference, opt$version)
HDF5Array::saveHDF5SummarizedExperiment(x = ref, dir = opt$output, verbose = opt$verbose,
										replace = opt$replace)
