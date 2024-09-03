#! /usr/bin/env Rscript


option_list <- list(
  optparse::make_option(c("-i", "--input"), type = "character",
              help="Input folder with 10X data"),
  optparse::make_option(c("-s", "--sample"), type = "character",
              help="Sample to subset")
)  

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

counts <- Seurat::Read10X(opt$input, gene.column = 1)
meta <- read.table(file.path(opt$input, "meta.tsv"), sep = "\t", header = TRUE)

for(sample in unique(meta$sample)) {
    cells <- meta[meta$sample == sample, ]$Cell_ID
    subset <- counts[, cells]
    DropletUtils::write10xCounts(file.path(getwd(), sample, "outs"), subset)
}
