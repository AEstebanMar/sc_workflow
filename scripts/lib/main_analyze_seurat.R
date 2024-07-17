##### Temporary docs from original functions
#' analyze_seurat
#' `analyze_seurat` performs all the analyses (individually or combining all
#' samples with and without integration)
#'
#' @param name sample name, or condition if integrate is TRUE
#' @param expermient experiment name
#' @param input directory with the single-cell data
#' @param output output directory (used when integrate is TRUE)
#' @param filter TRUE for using only detected cell-associated barcodes, FALSE
#' for using all detected barcodes
#' @param mincells min number of cells for which a feature is recorded
#' @param minfeats min number of features for which a cell is recorded
#' @param minqcfeats min number of features for which a cell is selected
#' @param percentmt max percentage of reads mapped to mitochondrial genes for 
#' which a cell is selected
#' @param normalmethod Normalization method
#' @param scalefactor Scale factor for cell-level normalization
#' @param hvgs Number of HVGs to be selected
#' @param ndims Number of PC to be used for clustering / UMAP / tSNE
#' @param resolution Granularity of the downstream clustering (higher values 
#' -> greater number of clusters)
#' @param integrate FALSE if we don't run integrative analysis, TRUE otherwise
#'
#' @keywords preprocessing, main
#' 
#' @return Seurat object

#' integrate_seurat
#' `integrate_seurat` allows subsetting a seurat object without requiring
#' literal strings.
#'
#' @param seu Merged eurat object
#' @param annotation_dir Directory with cluster annotation files
#' @param subset_column Column with categories to subset
#' @inheritParams analyze_seurat
#'
#' @returns A subsetted seurat object


main_analyze_seurat <- function(raw_seu, out_path = NULL, minqcfeats, percentmt,
                           		hvgs, ndims, resolution, dimreds_to_do,
                           		integrate = FALSE, cell_types_annotation,
                           		DEG_columns = NULL, int_columns = NULL){
  message('Normalizing data')
  merged_seu <- Seurat::NormalizeData(object = merged_seu, verbose = FALSE,
  									  normalization.method = normalmethod,
                         			  scale.factor = scalefactor)
  message('Scaling data')
  merged_seu <- Seurat::ScaleData(object = merged_seu, verbose = FALSE
  								  features = rownames(merged_seu))
  qc <- tag_qc(seu = raw_seu, minqcfeats = minqcfeats,percentmt = percentmt)
  seu <- subset(qc, subset = QC != 'High_MT,Low_nFeature')
  message('Finding variable features')
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = hvgs, verbose = FALSE)
  # dimreds_to_do: PCA/UMAP/tSNE
  message('Reducing dimensionality')  
  seu <- Seurat::RunPCA(seu, verbose = FALSE)
  if (integrate) {
  	message('Integrating seurat object')
  	seu <- harmony::RunHarmony(object = seu, "sample", plot_convergence = FALSE,
  	 							verbose = FALSE)
  	seu <- Seurat::RunUMAP(object = seu, dims = seq(1, ndims),
  						   reduction = "harmony", verbose = FALSE)
  	seu <- Seurat::FindNeighbors(object = seu, dims = seq(1, ndims),
  								 reduction = "harmony", verbose = FALSE)
  } else {
  	seu <- Seurat::RunUMAP(object = seu, dims = seq(1, ndims),
  						   reduction = "pca", verbose = FALSE)
  	seu <- Seurat::FindNeighbors(object = seu, dims = seq(1, ndims),
  								 reduction = "pca", verbose = FALSE)
  }
  if(length(int_columns) == 1) {
    markers <- get_sc_markers(seu = seu, cond = int_columns, DEG = FALSE, top = 200)
    markers <- collapse_markers(markers)
  }else{
    markers <- Seurat::FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  }
  if(clusters_annotation != "") {
  	message("Clusters annotation file provided. Annotating clusters")
  	seu <- annotate_clusters(seu = seu, new_clusters = clusters_annotation[[2]])
  }else if(cell_types_annotation != ""){
	message("Clusters annotation file not provided. Dynamically annotating
			   clusters.")
	anno_table <- read.table(opt$cell_types_annotation, sep = "\t", header = TRUE)
	annotated_clusters <- match_cell_types(markers, anno_table, top = 200)
	markers <- annotated_clusters$stats_table
	seu <- annotate_clusters(seu, annotated_clusters$cell_type)
  }else{
  	message("No cell type data provided. Clusters cannot be annotated")
  }
  message('Performing DEG analysis')
  seu <- SeuratObject::JoinLayers(seu)
  DEG_list <- list()
  if(DEG_columns == "") {
    DEG_conditions <- int_columns
  }else{
    DEG_conditions <- unlist(strsplit(DEG_columns, split = ","))
  }
  for(condition in DEG_conditions) {
    message(paste0("Calculating DEGs for condition ", condition))
    condition_DEGs <- get_sc_markers(seu = seu, cond = condition, DEG = TRUE)
    DEG_list[[condition]] <- condition_DEGs
  }
  if(!is.null(out_path)){
  	message('Writing results to disk')
  	saveRDS(seu, paste0(out_path, ".qc_seu.RDS"))
  	saveRDS(seu, paste0(out_path, ".seu.RDS"))
  	saveRDS(markers, paste0(out_path, ".markers.RDS")
  	saveRDS(DEG_list, paste0(out_path, ".DEGs.RDS")
  }
  final_results <- list()
  final_results$qc <- qc
  final_results$seu <- seu
  final_results$markers <- markers
  final_results$DEG_list <- DEG_list
  return(final_results)
}

