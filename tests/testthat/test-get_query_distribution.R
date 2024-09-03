test_that("get_query_distribution properly sums expression levels in samples", {
  library(Seurat)
  pbmc_smaller <- pbmc_small[, 1:15]
  pbmc_smaller@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
  query <- c("MS4A1", "CD79A", "HLA-DRB5")
  expected_pct <- data.frame(matrix(nrow = 3, ncol = 4))
  expected_pct[, 1] <- c("A", "B", "C")
  expected_pct[, 2] <- c(0, 0, 28.40)
  expected_pct[, 3] <- c(0, 0, 23.1)
  expected_pct[, 4] <- c(0, 5.06, 22.10)
  colnames(expected_pct) <- c("sample", query)
  output_pct <- get_query_distribution(pbmc_smaller, query, 3)
  expect_equal(output_pct, expected_pct)
})
