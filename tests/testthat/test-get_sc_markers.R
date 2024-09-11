library(Seurat)
test_that("get_sc_markers skips exclusive clusters in DEG analysis", {
  pbmc_smaller <- pbmc_small[, 1:15]
  pbmc_smaller@meta.data$groups <- "g1"
  pbmc_smaller@meta.data$groups[7:15] <- "g2"
  pbmc_smaller@meta.data$seurat_clusters <- 0
  pbmc_smaller@meta.data$seurat_clusters[3:5] <- 1
  pbmc_smaller@meta.data$seurat_clusters[15] <- 1
  expected_warning <- paste0("Cluster 2 contains less than three cells for ",
  							"condition 'g2'")
  expect_warning(suppressMessages(get_sc_markers(seu = pbmc_smaller,
  				 cond = "groups", DEG = TRUE, verbose = FALSE)), expected_warning)
  expect_false(suppressWarnings(
                suppressMessages(get_sc_markers(seu = pbmc_smaller,
  				      cond = "groups", DEG = TRUE, verbose = FALSE)[[2]][[1]])
                ))
})
