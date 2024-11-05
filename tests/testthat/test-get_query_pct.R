library(Seurat)
pbmc_smaller <- pbmc_small[, 1:15]
pbmc_smaller@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
pbmc_smaller@meta.data$seurat_clusters <- rep(0:4, 3)
cell_types <- c("0. typeA", "1. typeB", "2. typeC", "3. typeD", "4. typeE")
pbmc_smaller@meta.data$cell_type <- rep(cell_types, 3)
query <- c("MS4A1", "CD79A", "HLA-DRB5")

test_that("get_query_pct works in simple case", {
  query <- c("MS4A1", "CD79A", "HLA-DRB5")
  expected_df <- matrix(nrow = 3, ncol = 3)
  rownames(expected_df) <- c("A", "B", "C")
  expected_df[1, ] <- rep(0, 3)
  expected_df[2, ] <- c(0, 0, 20)
  expected_df[3, ] <- c(100, 80, 80)
  colnames(expected_df) <- query
  output_df <- get_query_pct(pbmc_smaller, query, "sample")
  expect_equal(output_df, expected_df)
})

test_that("get_query_pct gives warning if any query genes are not found", {
  missing_query <- c("MS4A1", "CD79A", "HLA-DRB5", "NOEXPA", "NOEXPB")
  expected_df <- matrix(nrow = 3, ncol = 3)
  rownames(expected_df) <- c("A", "B", "C")
  expected_df[1, ] <- rep(0, 3)
  expected_df[2, ] <- c(0, 0, 20)
  expected_df[3, ] <- c(100, 80, 80)
  colnames(expected_df) <- query
  expect_warning(get_query_pct(pbmc_smaller, missing_query, "sample"), "NOEXPA, NOEXPB")
  output_df <- suppressWarnings(get_query_pct(pbmc_smaller, missing_query, "sample"))
  expect_equal(output_df, expected_df)
})

test_that("get_query_pct works with query of length one", {
  single_query <- c("HLA-DRB5")
  expected_df <- matrix(nrow = 3, ncol = 1)
  rownames(expected_df) <- c("A", "B", "C")
  expected_df[, 1] <- c(0, 20, 80)
  colnames(expected_df) <- single_query
  output_df <- get_query_pct(pbmc_smaller, single_query, "sample")
  expect_equal(output_df, expected_df)
})

test_that("get_query_pct works with alternate 'by' arguments", {
  single_query <- c("HLA-DRB5")
  expected_df <- matrix(nrow = 5, ncol = 1)
  rownames(expected_df) <- c(paste0(0:4, ". type", toupper(letters[1:5])))
  expected_df[, 1] <- c(33, 67, 33, 0, 33)
  colnames(expected_df) <- single_query
  output_df <- get_query_pct(pbmc_smaller, single_query, "cell_type")
  testthat::expect_equal(output_df, expected_df)
})

# Test disabled until I figure out a way to add a "counts" element to assays
# in a way that lets me access it in the same way that get_query_pct does.
# That method breaks in test dataset, I have not figured out what makes it
# different. Does not have to do with different seurat version.
# test_that("get_query_pct works with 'by' argument of length 2", {
#   pbmc_updated <- Seurat::CreateSeuratObject(counts = pbmc_smaller$RNA$counts)
#   pbmc_updated@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
#   pbmc_updated@meta.data$seurat_clusters <- rep(0:4, 3)
#   types <- paste0("type", toupper(letters[1:5]))
#   matA <- matrix(data = 0, nrow = 5, ncol = 3)
#   matB = matC <- matA
#   matB[, 3] <- c(0, 100, rep(0, 3))
#   matC[, 1] <- rep(100, 5)
#   matC[, 2] <- c(0, rep(100, 4))
#   matC[, 3] <- c(rep(100, 3), 0, 100)
#   colnames(matA) = colnames(matB) = colnames(matC) <- query
#   rownames(matA) = rownames(matB) = rownames(matC) <- types
#   expected_list <- list(A = matA, B = matB, C = matC)
#   output_list <- get_query_pct(pbmc_updated, query,
#                                by = c("sample", "seurat_clusters"))
#   expect_equal(output_list, expected_list)
# })
