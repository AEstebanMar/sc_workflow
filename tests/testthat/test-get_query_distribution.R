
test_that("get_query_distribution properly sums expression levels in samples", {
  pbmc_tiny <- pbmc_small[, 1:15]
  pbmc_tiny@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
  query <- c("MS4A1", "CD79A", "HLA-DRB5")
  expected_exp <- data.frame(matrix(nrow = 3, ncol = 3))
  rownames(expected_exp) <- c("A", "B", "C")
  expected_exp[, 1] <- c(0, 0, 28.40)
  expected_exp[, 2] <- c(0, 0, 23.1)
  expected_exp[, 3] <- c(0, 5.06, 22.10)
  colnames(expected_exp) <- query
  output_exp <- get_query_distribution(pbmc_tiny, query, 3)
  expect_equal(output_exp, expected_exp)
})
