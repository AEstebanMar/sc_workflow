test_that("Cell type matching sorted simple case"),{
	markers_df <- data.frame(samples = seq(1, 4),
							 cluster = seq(1, 4))
	rownames(markers_df) <- c("A", "B", "C", "D")
	anno_table <- data.frame(markers = c("A", "B", "C"),
							 type = c("Type1", "Type2", "Type3"))
	expected_df <- markers_df
	expected_df$cell_type <- c("Type1", "Type2", "Type3", "Unknown")
	expected_df$gene <- rownames(expected_df)
	expected_df <- expected_df[order(expected_df$cluster), ]
	output_df <- match_cell_types(markers_df, anno_table)
	expect_equal(output_df$stats_table, expected_df)
}

test_that("Cell type matching unsorted simple case"),{
	markers_df <- data.frame(samples = c(4, 3, 2, 1),
							 cluster = c(2, 3, 1, 4))
	rownames(markers_df) <- c("A", "C", "B", "D")
	anno_table <- data.frame(markers = c("C", "B", "A"),
							 type = c("Type3", "Type2", "Type1"))
	expected_df <- markers_df
	expected_df$cell_type <- c("Type1", "Type3", "Type2", "Unknown")
	expected_df$gene <- rownames(expected_df)
	expected_df <- expected_df[order(expected_df$cluster), ]
	output_df <- match_cell_types(markers_df, anno_table)
	expect_equal(output_df$stats_table, expected_df)
}

test_that("Cell type matching complex case",{
	markers_df <- data.frame(samples = seq(1, 50),
							 cluster = c(rep(1, 20), rep(2, 20), rep (3, 10)))
	genes <- paste0("gene", seq(1, 50))
	rownames(markers_df) <- genes
	types <- c(rep("type1", 20), rep("type2", 13), rep("type3", 17))
	anno_table <- data.frame(markers = genes,
							 type = types)
	expected_df <- markers_df
	expected_df$cell_type <-  c(rep("type1", 20), rep("type3", 20), rep("type2", 10))
	expected_df$gene <- rownames(expected_df)
	expected_df <- expected_df[order(expected_df$cluster), ]
	output_df <- match_cell_types(markers_df, anno_table)
	expect_equal(output_df$stats_table, expected_df)
})