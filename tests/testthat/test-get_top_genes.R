library(Seurat)
pbmc_smaller <- pbmc_small[, 1:15]
pbmc_smaller@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
pbmc_smaller@meta.data$seurat_clusters <- rep(0:4, 3)
cell_types <- c("0. typeA", "1. typeB", "2. typeC", "3. typeD", "4. typeE")
pbmc_smaller@meta.data$cell_type <- rep(cell_types, 3)
query <- c("MS4A1", "CD79A", "HLA-DRB5")
counts <- GetAssayData(pbmc_smaller)
gene_list <- c("MS4A1", "CD79B", "CD79A", "LTB", "TCL1A")
expected_vector <- c("MS4A1", "CD79B", "CD79A","LTB", "TCL1A")
counts <- counts[gene_list, ]
pbmc_mini <- subset(pbmc_smaller, features = rownames(counts))

test_that("get_top_genes works as expected", {
	simple_vector <- get_top_genes(pbmc_mini, top = 5)
	greater_vector <- get_top_genes(pbmc_mini, top = 10e99)
	expect_equal(expected_vector, simple_vector)
	expect_equal(expected_vector, greater_vector)
	expect_error(get_top_genes(pbmc_mini, top = 0), "greater than 1, was 0")
})

test_that("get_top_genes works even if N is greater than number of expressed
		   genes", {
	output_vector <- get_top_genes(pbmc_mini, top = 10e99999)
	testthat::expect_equal(expected_vector, output_vector)
})
