test_that("get_query_pct works in simple case", {
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
  output_df <- get_query_pct(pbmc_smaller, query)
  expect_equal(output_df, expected_df)
})

test_that("get_query_pct gives warning if any query genes are not found", {
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
  expect_warning(get_query_pct(pbmc_smaller, query), "NOEXPA, NOEXPB")
  output_df <- suppressWarnings(get_query_pct(pbmc_smaller, query))
  expect_equal(output_df, expected_df)
})

test_that("get_query_pct works with query of length one", {
  library(Seurat)
  pbmc_smaller <- pbmc_small[, 1:15]
  pbmc_smaller@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
  query <- c("HLA-DRB5")
  expected_df <- matrix(nrow = 3, ncol = 1)
  rownames(expected_df) <- c("A", "B", "C")
  expected_df[, 1] <- c(0, 20, 80)
  colnames(expected_df) <- query
  output_df <- get_query_pct(pbmc_smaller, query)
  expect_equal(output_df, expected_df)
})
