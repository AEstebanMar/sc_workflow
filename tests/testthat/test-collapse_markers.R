
test_that("collapse_markers base use case", {
  markers_1 <- data.frame(a_avg_log2FC = rep(1, 3), b_avg_log2FC = rep (-1, 3))
  markers_2 <- data.frame(a_avg_log2FC = rep(-2, 3), b_avg_log2FC = rep (2, 3))
  markers_list <- list(cluster_1 = markers_1, cluster_2 = markers_2)
  expected_df <-
  output_df <- collapse_markers(markers_list)
  expected_df <- data.frame(gene = as.character(rep(1:3, 2)),
                            a_avg_log2FC = c(rep(1, 3), rep(-2, 3)),
                            b_avg_log2FC = c(rep(-1, 3), rep (2, 3)),
                            cluster = c(rep(0, 3), rep(1, 3)))
  expect_equal(output_df, expected_df)
})
