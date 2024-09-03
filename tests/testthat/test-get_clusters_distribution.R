test_that("get_clusters_distribution properly calculates percentages", {
  library(Seurat)
  pbmc_smaller <- pbmc_small[, 1:15]
  pbmc_smaller@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
  pbmc_smaller@meta.data$seurat_clusters <- c(rep(1, 3), rep(2, 10), rep(3, 2))
  expected_pct <- matrix(nrow = 3, ncol = 3)
  expected_pct[1, ] <- c(0.6, 0.4, 0)
  expected_pct[2, ] <- c(0, 1.0, 0)
  expected_pct[3, ] <- c(0, 0.6, 0.4)
  expected_pct <- expected_pct * 100
  output_pct <- get_clusters_distribution(pbmc_smaller, 2)
  dimnames(output_pct) <- NULL
  output_pct <- apply(output_pct, 2, rbind)
  expect_equal(output_pct, expected_pct)
})