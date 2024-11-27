
gen_pbmc_tiny <- pbmc_small[, 1:15]
gen_pbmc_tiny@meta.data$groups <- "g1"
gen_pbmc_tiny@meta.data$groups[7:15] <- "g2"
gen_pbmc_tiny@meta.data$seurat_clusters <- 0
gen_pbmc_tiny@meta.data$seurat_clusters[3:5] <- 1
gen_pbmc_tiny@meta.data$seurat_clusters[15] <- 1
gen_pbmc_tiny@meta.data$cell_types <- gen_pbmc_tiny@meta.data$seurat_clusters

test_that("get_sc_markers skips exclusive clusters in DEG analysis", {
  expected_warning <- paste0("Cluster 2 contains less than three cells for ",
  							"condition 'g2'")
  expect_warning(suppressMessages(get_sc_markers(seu = gen_pbmc_tiny,
  				       cond = "groups", DEG = TRUE, verbose = FALSE,
                 subset_by = "seurat_clusters")), expected_warning)
  expect_false(suppressWarnings(
                suppressMessages(get_sc_markers(seu = gen_pbmc_tiny,
  				      cond = "groups", DEG = TRUE, verbose = FALSE,
                subset_by = "seurat_clusters")$markers[[2]][[1]])
                ))
})

test_that("get_sc_markers skips exclusive clusters in DEG analysis, alternate
           idents", {
  expected_warning <- paste0("Cluster 2 contains less than three cells for ",
                "condition 'g2'")
  expect_warning(suppressMessages(get_sc_markers(seu = gen_pbmc_tiny,
                 cond = "groups", DEG = TRUE, verbose = FALSE,
                 subset_by = "cell_types")), expected_warning)
  expect_false(suppressWarnings(
                suppressMessages(get_sc_markers(seu = gen_pbmc_tiny,
                cond = "groups", DEG = TRUE, verbose = FALSE,
                subset_by = "cell_types")$markers[[2]][[1]])
                ))
})
