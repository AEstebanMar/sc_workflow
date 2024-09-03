test_that("Cell type matching sorted simple case",{
	markers_df <- data.frame(samples = seq(1, 4),
							 cluster = seq(1, 4),
							 gene = c("A", "B", "C", "D"))
	markers_df$p_val_adj <- rep(1e-5, nrow(markers_df))
	markers_df$avg_log2FC <- rep(1, nrow(markers_df))
	rownames(markers_df) <- c("A", "B", "C", "D")
	cell_annotation <- data.frame(markers = c("A", "B", "C"),
							 type = c("Type1", "Type2", "Type3"))
	expected_df <- markers_df
	expected_df$cell_type <- c("1. Type1", "2. Type2", "3. Type3", "4. Unknown")
	expected_df$gene <- rownames(expected_df)
	expected_df <- expected_df[order(expected_df$cluster), ]
	output_df <- match_cell_types(markers_df = markers_df,
								  cell_annotation = cell_annotation)$stats_table
	expect_equal(output_df, expected_df)
})

test_that("Cell type matching unsorted simple case",{
	markers_df <- data.frame(samples = c(4, 3, 2, 1),
							 cluster = c(2, 3, 1, 4),
							 gene = c("A", "C", "B", "D"))
	markers_df$p_val_adj <- rep(1e-5, nrow(markers_df))
	markers_df$avg_log2FC <- rep(1, nrow(markers_df))
	rownames(markers_df) <- c("A", "C", "B", "D")
	cell_annotation <- data.frame(markers = c("C", "B", "A"),
							 type = c("Type3", "Type2", "Type1"))
	expected_df <- markers_df
	expected_df$cell_type <- c("2. Type1", "3. Type3", "1. Type2", "4. Unknown")
	expected_df$gene <- rownames(expected_df)
	expected_df <- expected_df[order(expected_df$cluster), ]
	output_df <- match_cell_types(markers_df, cell_annotation)$stats_table
	expect_equal(output_df, expected_df)
})

test_that("Cell type matching complex case",{
	markers_df <- data.frame(samples = seq(1, 50),
							 cluster = c(rep(1, 20), rep(2, 20), rep (3, 10)))
	markers_df$p_val_adj <- rep(1e-5, nrow(markers_df))
	markers_df$avg_log2FC <- rep(1, nrow(markers_df))
	genes <- paste0("gene", seq(1, 50))
	rownames(markers_df) <- genes
	markers_df$gene <- genes
	types <- c(rep("type1", 20), rep("type2", 13), rep("type3", 17))
	cell_annotation <- data.frame(markers = genes,
							 type = types)
	expected_df <- markers_df
	expected_df$cell_type <-  c(rep("1. type1", 20), rep("2. type2", 20),
								rep("3. type3", 10))
	expected_df$gene <- rownames(expected_df)
	expected_df <- expected_df[order(expected_df$cluster), ]
	output_df <- match_cell_types(markers_df, cell_annotation)$stats_table
	expect_equal(output_df, expected_df)
})