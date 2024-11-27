
gen_pbmc_tiny <- pbmc_small[, 1:15]
gen_pbmc_tiny@meta.data$groups <- "g1"
gen_pbmc_tiny@meta.data$groups[7:15] <- "g2"
gen_pbmc_tiny@meta.data$seurat_clusters <- 0
gen_pbmc_tiny@meta.data$seurat_clusters[1] <- 1
gen_pbmc_tiny@meta.data$seurat_clusters[15] <- 1
gen_pbmc_tiny@meta.data$cell_types <- "typeA"
gen_pbmc_tiny@meta.data$cell_types[3:5] <- "typeB"
gen_pbmc_tiny@meta.data$cell_types[15] <- "typeB"
test_that("has_exclusive_idents works in base case", {
  expected_warning <- "Defaulting to general marker analysis"
  expect_warning(has_exclusive_idents(seu = gen_pbmc_tiny, cond = "groups",
                                  idents = "seurat_clusters"), expected_warning)
  expect_true(suppressWarnings(has_exclusive_idents(seu = gen_pbmc_tiny,
                                  idents = "seurat_clusters", cond = "groups")))
})

test_that("has_exclusive_idents works with alternate ident", {
  expected_warning <- "Defaulting to general marker analysis"
  expect_warning(has_exclusive_idents(seu = gen_pbmc_tiny, cond = "groups",
                                  idents = "cell_types"), expected_warning)
  expect_true(suppressWarnings(has_exclusive_idents(seu = gen_pbmc_tiny,
                                  idents = "cell_types", cond = "groups")))
})

test_that("has_exclusive_idents can handle more than one positive", {
  expected_warning <- paste0("g1-1, g2-1")
  expect_warning(has_exclusive_idents(seu = gen_pbmc_tiny, cond = "groups",
                                  idents = "seurat_clusters"), expected_warning)
  expect_true(suppressWarnings(has_exclusive_idents(seu = gen_pbmc_tiny,
                                  idents = "seurat_clusters", cond = "groups")))
})

test_that("has_exclusive_idents can predict pairs that should appear but
           do not (behaves correctly for strictly exclusive idents)", {
  pbmc_tiny <- pbmc_small[, 1:15]
  pbmc_tiny@meta.data$groups <- "g1"
  pbmc_tiny@meta.data$groups[7:15] <- "g2"
  pbmc_tiny@meta.data$seurat_clusters <- 0
  pbmc_tiny@meta.data$seurat_clusters[7:15] <- 1
  expected_warning <- paste0("g2-0, g1-1")
  expect_warning(has_exclusive_idents(seu = pbmc_tiny, cond = "groups",
                                  idents = "seurat_clusters"), expected_warning)
  expect_true(suppressWarnings(has_exclusive_idents(seu = pbmc_tiny,
                                  idents = "seurat_clusters", cond = "groups")))
})

test_that("has_exclusive_idents can identify that no exclusive idents
           exist", {
  pbmc_tiny <- pbmc_small[, 1:15]
  pbmc_tiny@meta.data$seurat_clusters <- 0
  pbmc_tiny@meta.data$seurat_clusters[c(8:15)] <- 1
  expect_false(suppressWarnings(has_exclusive_idents(seu = pbmc_tiny,
                                  idents = "seurat_clusters", cond = "groups")))
})
