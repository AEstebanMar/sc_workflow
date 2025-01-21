#! /usr/bin/env Rscript
# Custom test file for non-package structure


library(testthat)
library(Seurat)
devtools::load_all("~/dev_R/ExpHunterSuite")
test_dir("./testthat/")
