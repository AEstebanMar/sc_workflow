test_that("annotate_clusters properly marks duplicated cell types", {
  library(Seurat)
  pbmc_smaller <- pbmc_small[, 1:15]
  Idents(pbmc_smaller) <- 1:15
  new_clusters <- rep("TypeA", 15)
  new_clusters[c(3, 7, 10)] <- "TypeB"
  new_clusters[c(2, 5, 13)] <- "TypeC"
  annotated <- suppressWarnings(annotate_clusters(pbmc_smaller, new_clusters))
  output <- annotated@meta.data$named_clusters
  expected <- c("TypeA (a)", "TypeC (a)", "TypeB (a)", "TypeA (b)", "TypeC (b)",
                "TypeA (c)", "TypeB (b)", "TypeA (d)", "TypeA (e)", "TypeB (c)",
                "TypeA (f)", "TypeA (g)", "TypeC (c)", "TypeA (h)", "TypeA (i)")
  expect_equal(as.character(output), expected)
})
