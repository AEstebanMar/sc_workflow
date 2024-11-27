
test_that("Cell type matching sorted simple case", {
	markers_df <- data.frame(samples = seq(1, 4),
							 cluster = seq(1, 4),
							 gene = c("A", "B", "C", "D"))
	markers_df$p_val_adj <- rep(1e-5, nrow(markers_df))
	markers_df$avg_log2FC <- rep(1, nrow(markers_df))
	cell_annotation <- data.frame(markers = c("A", "B", "C"),
							 type = c("Type1", "Type2", "Type3"))
	expected_df <- markers_df
	expected_df$cell_type <- c("1. Type1", "2. Type2", "3. Type3", "4. Unknown")
	expected_df <- expected_df[order(expected_df$cluster), ]
	output_df <- match_cell_types(markers_df = markers_df,
								  cell_annotation = cell_annotation)$stats_table
	expect_equal(output_df, expected_df)
})

test_that("Cell type matching unsorted simple case", {
	markers_df <- data.frame(samples = c(4, 3, 2, 1),
							 cluster = c(2, 3, 1, 4),
							 gene = c("A", "C", "B", "D"))
	markers_df$p_val_adj <- rep(1e-5, nrow(markers_df))
	markers_df$avg_log2FC <- rep(1, nrow(markers_df))
	cell_annotation <- data.frame(markers = c("C", "B", "A"),
							 type = c("Type3", "Type2", "Type1"))
	expected_df <- markers_df
	expected_df$cell_type <- c("2. Type1", "3. Type3", "1. Type2", "4. Unknown")
	expected_df <- expected_df[order(expected_df$cluster), ]
	output_df <- match_cell_types(markers_df, cell_annotation)$stats_table
	expect_equal(output_df, expected_df)
})

test_that("Cell type matching complex case", {
	markers_df <- data.frame(samples = seq(1, 50),
							 cluster = c(rep(1, 20), rep(2, 20), rep (3, 10)))
	markers_df$p_val_adj <- rep(1e-5, nrow(markers_df))
	markers_df$avg_log2FC <- rep(1, nrow(markers_df))
	genes <- paste0("gene", seq(1, 50))
	markers_df$gene <- genes
	types <- c(rep("type1", 20), rep("type2", 13), rep("type3", 17))
	cell_annotation <- data.frame(markers = genes, type = types)
	expected_df <- markers_df
	expected_df$cell_type <-  c(rep("1. type1", 20), rep("2. type2", 20),
								rep("3. type3", 10))
	expected_df <- expected_df[order(expected_df$cluster), ]
	output_df <- match_cell_types(markers_df, cell_annotation)$stats_table
	expect_equal(output_df, expected_df)
})

test_that("More than one cluster assigned the same cell type", {
	markers_df <- data.frame(samples = 1:4, cluster = 1:4)
	markers_df$p_val_adj <- rep(1e-5, 4)
	markers_df$avg_log2FC <- rep(1, 4)
	genes <- paste0("gene", seq(3))
	genes <- c(genes, paste0("gene", 1))
	markers_df$gene <- genes
	types <- c("type1", "type2", "type3")
	cell_annotation <- data.frame(markers = genes[1:12], type = types)
	expected_df <- markers_df
	expected_df$cell_type <- c("1. type1 (a)", "2. type2", "3. type3",
							   "4. type1 (b)")
	output_df <- match_cell_types(markers_df, cell_annotation)$stats_table
	expect_equal(output_df, expected_df)
})

test_that("Draw between more than one cell type", {
	markers_df <- data.frame(samples = 1:4, cluster = 1:4)
	markers_df$p_val_adj <- rep(1e-5, 4)
	markers_df$avg_log2FC <- rep(1, 4)
	genes <- paste0("gene", c(1, 1, 2))
	genes <- c(genes, paste0("gene", 1))
	markers_df$gene <- genes
	types <- c("type1", "type2", "type3")
	cell_annotation <- data.frame(markers = genes[1:12], type = types)
	expected_df <- markers_df
	expected_df$cell_type <- c("1. type1 / type2 (a)", "2. type1 / type2 (b)",
							   "3. type3", "4. type1 / type2 (c)")
	output_df <- match_cell_types(markers_df, cell_annotation)$stats_table
	expect_equal(output_df, expected_df)
})
