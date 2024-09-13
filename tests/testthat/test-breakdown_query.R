
test_that("breakdown_query works in simple case", {
  library(Seurat)
  pbmc_smaller <- pbmc_small[, 1:15]
  pbmc_smaller@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
  query <- c("MS4A1", "CD79A", "HLA-DRB5")
  expected_df <- matrix(nrow = 3, ncol = 3)
  rownames(expected_df) <- c("A", "B", "C")
  expected_df[1, ] <- rep(0, 3)
  expected_df[2, ] <- c(0, 0, 20)
  expected_df[3, ] <- c(100, 80, 80)
  colnames(expected_df) <- query
  output_df <- breakdown_query(pbmc_smaller, query, "sample")
  expect_equal(output_df, expected_df)
})

test_that("breakdown_query gives warning if any query genes are not found", {
  library(Seurat)
  pbmc_smaller <- pbmc_small[, 1:15]
  pbmc_smaller@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
  query <- c("MS4A1", "CD79A", "HLA-DRB5", "NOEXPA", "NOEXPB")
  expected_df <- matrix(nrow = 3, ncol = 3)
  rownames(expected_df) <- c("A", "B", "C")
  expected_df[1, ] <- rep(0, 3)
  expected_df[2, ] <- c(0, 0, 20)
  expected_df[3, ] <- c(100, 80, 80)
  colnames(expected_df) <- query[1:3]
  expect_warning(breakdown_query(pbmc_smaller, query, "sample"), "NOEXPA, NOEXPB")
  output_df <- suppressWarnings(breakdown_query(pbmc_smaller, query, "sample"))
  expect_equal(output_df, expected_df)
})

test_that("breakdown_query works with query of length one", {
  library(Seurat)
  pbmc_smaller <- pbmc_small[, 1:15]
  pbmc_smaller@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
  query <- c("HLA-DRB5")
  expected_df <- matrix(nrow = 3, ncol = 1)
  rownames(expected_df) <- c("A", "B", "C")
  expected_df[, 1] <- c(0, 20, 80)
  colnames(expected_df) <- query
  output_df <- breakdown_query(pbmc_smaller, query, "sample")
  expect_equal(output_df, expected_df)
})

test_that("breakdown_query works with alternate 'by' arguments", {
  library(Seurat)
  pbmc_smaller <- pbmc_small[, 1:15]
  pbmc_smaller@meta.data$seurat_clusters <- c(rep(0, 5), rep (1, 5), rep (2, 5))
  query <- c("HLA-DRB5")
  expected_df <- matrix(nrow = 3, ncol = 1)
  rownames(expected_df) <- c(0:2)
  expected_df[, 1] <- c(0, 20, 80)
  colnames(expected_df) <- query
  output_df <- breakdown_query(pbmc_smaller, query, "seurat_clusters")
  testthat::expect_equal(output_df, expected_df)
})
