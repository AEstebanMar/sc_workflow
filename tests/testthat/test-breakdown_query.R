
library(Seurat)
pbmc_smaller <- pbmc_small[, 1:15]
pbmc_smaller@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
pbmc_smaller@meta.data$seurat_clusters <- c(rep(0, 5), rep (1, 5), rep (2, 5))
query <- c("MS4A1", "CD79A", "HLA-DRB5")
test_that("breakdown_query works in simple case", {
  expected_df <- c(0.33, 0.27, 0.33)
  names(expected_df) <- query
  output_df <- signif(breakdown_query(pbmc_smaller, query), 2)
  expect_equal(output_df, expected_df)
})

test_that("breakdown_query gives warning if any query genes are not found", {
  library(Seurat)
  missing_query <- c(query, "NOEXPA", "NOEXPB")
  expected_df <- c(0.33, 0.27, 0.33)
  names(expected_df) <- query
  expect_warning(breakdown_query(pbmc_smaller, missing_query), "NOEXPA, NOEXPB")
  output_df <- suppressWarnings(signif(breakdown_query(pbmc_smaller,
                                                       missing_query), 2))
  expect_equal(output_df, expected_df)
})

test_that("breakdown_query works with query of length one", {
  single_query <- c("HLA-DRB5")
  expected_df <- 0.33
  names(expected_df) <- single_query
  output_df <- signif(breakdown_query(pbmc_smaller, single_query), 2)
  expect_equal(output_df, expected_df)
})
