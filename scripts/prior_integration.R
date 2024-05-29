#! /usr/bin/env Rscript


##########################################
## LOAD LIBRARIES
##########################################
# Obtain this script directory
full.fpath <- normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                commandArgs())], '='))[2])

main_path_script <- dirname(full.fpath)
root_path <- file.path(main_path_script)
# Load custom libraries
# devtools::load_all(file.path(root_path))

source_folder <- file.path(root_path, 'lib')
library(Seurat)
library(scCustomize)
source(file.path(source_folder, "preprocessing_library.R"))

##########################################
## OPTPARSE
##########################################

option_list <- list(
  optparse::make_option(c("-d", "--exp_design"), type = "character",
              help="Input file with the experiment design"),
  optparse::make_option(c("-o", "--output"), type = "character",
              help="Output folder"),
  optparse::make_option(c("-c", "--condition"), type = "character",
              help="Name of the column in the experimental design file we want to use for integration"),
  optparse::make_option(c("-i", "--integration_file"), type = "character",
              help="Name of file that will contain integration subset names"),
  optparse::make_option(c("-e", "--experiment_name"), type = "character",
              help="Experiment name"),
  optparse::make_option(c("--count_path"), type = "character",
              help="Count results folder"),
  optparse::make_option(c("--suffix"), type = "character",
              help="Suffix to specific file")
)  


opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

############
### Main ###
############

# Divide sample names by the condition

exp_subsets <- do_subsetting(exp_design = opt$exp_design,
	                           column = opt$condition)


# Add them to a file so we pass it to the SAMPLES_FILE variable in autoflow_launcher.sh

if (file.exists(opt$integration_file)) { # Integration subsets file creation (must be empty if already exists from previous runs)
  file_conn <- file(opt$integration_file, open = "w")
  close(file_conn)
} else {
  file.create(opt$integration_file)
}

cat(names(exp_subsets), file = opt$integration_file, sep = "\n", append = FALSE) # the second quotes are for having an extra empty line at the end (see README in Github for details)

# Create and save the "before" (with no preprocessing) Seurat object

for (cond in names(exp_subsets)){
  folder_name <- opt$output
  if (!file.exists(folder_name)) {
    dir.create(folder_name)
  }
  seu <- merge_condition(exp_cond = cond,
                         samples = exp_subsets[[cond]],
                         exp_design = opt$exp_design,
                         count_path = opt$count_path,
                         suffix = opt$suffix)
  saveRDS(seu, file = file.path(folder_name, paste0(opt$experiment_name, ".", cond, ".before.seu.RDS")))
}