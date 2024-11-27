
test_that("annotate_clusters simply assigns names to clusters", {
  pbmc_tiny <- pbmc_small[, 1:15]
  Idents(pbmc_tiny) <- 1:15
  new_clusters <- rep("TypeA", 15)
  new_clusters[c(3, 7, 10)] <- "TypeB"
  new_clusters[c(2, 5, 13)] <- "TypeC"
  annotated <- suppressWarnings(annotate_clusters(pbmc_tiny, new_clusters))
  output <- annotated@meta.data$cell_type
  expected <- c("TypeA", "TypeC", "TypeB", "TypeA", "TypeC",
                "TypeA", "TypeB", "TypeA", "TypeA", "TypeB",
                "TypeA", "TypeA", "TypeC", "TypeA", "TypeA")
  expect_equal(as.character(output), expected)
})
