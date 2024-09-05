library(Seurat)
test_that("has_exclusive_clusters works in base case"), {
  
  pbmc_smaller <- pbmc_small[, 1:15]
  pbmc_smaller@meta.data$groups <- "g1"
  pbmc_smaller@meta.data$groups[7:15] <- "g2"
  pbmc_smaller@meta.data$seurat_clusters <- 0
  pbmc_smaller@meta.data$seurat_clusters[3:5] <- 1
  pbmc_smaller@meta.data$seurat_clusters[15] <- 1
  expected_warning <- "Defaulting to general marker analysis"
  expect_warning(has_exclusive_clusters(seu = pbmc_smaller, cond = "groups"),
                 expected_warning)
  expect_true(suppressWarnings(has_exclusive_clusters(seu = pbmc_smaller,
                                                      cond = "groups")))
}

test_that("has_exclusive_clusters can handle more than one positive"), {
  pbmc_smaller <- pbmc_small[, 1:15]
  pbmc_smaller@meta.data$groups <- "g1"
  pbmc_smaller@meta.data$groups[7:15] <- "g2"
  pbmc_smaller@meta.data$seurat_clusters <- 0
  pbmc_smaller@meta.data$seurat_clusters[1] <- 1
  pbmc_smaller@meta.data$seurat_clusters[15] <- 1
  expected_warning <- paste0("g1-1, g2-1")
  expect_warning(has_exclusive_clusters(seu = pbmc_smaller, cond = "groups"),
                 expected_warning)
  expect_true(suppressWarnings(has_exclusive_clusters(seu = pbmc_smaller,
                                                      cond = "groups")))
}

test_that("has_exclusive_clusters can predict pairs that should appear but
           do not (behaves correctly for strictly exclusive clusters)"), {
  pbmc_smaller <- pbmc_small[, 1:15]
  pbmc_smaller@meta.data$groups <- "g1"
  pbmc_smaller@meta.data$groups[7:15] <- "g2"
  pbmc_smaller@meta.data$seurat_clusters <- 0
  pbmc_smaller@meta.data$seurat_clusters[7:15] <- 1
  expected_warning <- paste0("g2-0, g1-1")
  expect_warning(has_exclusive_clusters(seu = pbmc_smaller, cond = "groups"),
                 expected_warning)
  expect_true(suppressWarnings(has_exclusive_clusters(seu = pbmc_smaller,
                                                      cond = "groups")))
}

test_that("has_exclusive_clusters can identify that no exclusive clusters
           exist"), {
  pbmc_smaller <- pbmc_small[, 1:15]
  pbmc_smaller@meta.data$seurat_clusters <- 0
  pbmc_smaller@meta.data$seurat_clusters[c(8:15)] <- 1
  expect_false(suppressWarnings(has_exclusive_clusters(seu = pbmc_smaller,
                                                       cond = "groups")))
}
