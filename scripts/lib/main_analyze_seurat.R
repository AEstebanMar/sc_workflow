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
#' @param seu Merged seurat object
#' @param annotation_dir Directory with cluster annotation files
#' @param subset_column Column with categories to subset
#' @param sigfig Significant decimal figures cutoff in query and cluster
#' distribution analysis.
#' @param query Vector of query genes whose expression to analyse.
#' @returns A subsetted seurat object

main_analyze_seurat <- function(seu, minqcfeats, percentmt, query, sigfig = 2,
                           		  resolution, dimreds_to_do, p_adj_cutoff = 5e-3,
                           		  integrate = FALSE, cluster_annotation = NULL,
                           		  cell_annotation, DEG_columns = NULL,
                           		  scalefactor = 10000, hvgs, int_columns = NULL,
                           		  normalmethod = "LogNormalize", ndims,
                                verbose = FALSE, output = getwd(),
                                save_RDS = FALSE, reduce = FALSE){
  colnames(seu@meta.data) <- tolower(colnames(seu@meta.data))
  qc <- tag_qc(seu = seu, minqcfeats = minqcfeats, percentmt = percentmt)
  if(!reduce) {
    seu <- subset(qc, subset = qc != 'High_MT,Low_nFeature')
  } else {
    message("Reduce argument is set to TRUE. Skipping QC subsetting")
    seu <- qc
  }
  message('Normalizing data')
  seu <- Seurat::NormalizeData(object = seu, verbose = verbose,
  									           normalization.method = normalmethod,
                               scale.factor = scalefactor)
  message('Scaling data')
  seu <- Seurat::ScaleData(object = seu, verbose = verbose,
  								         features = rownames(seu))
  message('Finding variable features')
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = hvgs, verbose = verbose)
  message('Reducing dimensionality')  
  seu <- Seurat::RunPCA(seu, verbose = verbose)
  reduction <- "pca"
  if(integrate) {
  	message('Integrating seurat object')
  	seu <- harmony::RunHarmony(object = seu, "sample", plot_convergence = FALSE,
  	 							verbose = verbose)
    reduction <- "harmony"
  }
  seu <- Seurat::RunUMAP(object = seu, dims = seq(1, ndims),
                 reduction = reduction, verbose = verbose)
  seu <- Seurat::FindNeighbors(object = seu, dims = seq(1, ndims),
                   reduction = reduction, verbose = verbose)
  seu <- Seurat::FindClusters(seu, resolution = resolution, verbose = verbose)
  seu <- SeuratObject::JoinLayers(seu)
  run_conserved <- ifelse(test = length(int_columns) == 1, no = FALSE,
                          yes = !has_exclusive_clusters(seu = seu,
                                                        cond = int_columns))
  if(run_conserved) {
    markers <- get_sc_markers(seu = seu, cond = int_columns, DEG = FALSE)
    markers <- collapse_markers(markers)
  }else{
    markers <- Seurat::FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25,
                                      logfc.threshold = 0.25, verbose = verbose)
    rownames(markers) <- NULL
  }
  if(!is.null(cluster_annotation)) {
  	message("Clusters annotation file provided. Annotating clusters")
  	seu <- annotate_clusters(seu = seu, new_clusters = cluster_annotation$name)
  }else if(!is.null(cell_annotation)){
	  message("Clusters annotation file not provided. Dynamically annotating
			       clusters.")
	  annotated_clusters <- match_cell_types(markers_df = markers,
                                           cell_annotation = cell_annotation,
                                           p_adj_cutoff = p_adj_cutoff)
	  markers <- annotated_clusters$summary
	  seu <- annotate_clusters(seu, annotated_clusters$cell_types)
  }else{
  	warning("No data provided for cluster annotation")
  }
  clusters_pct <- get_clusters_distribution(seu = seu, sigfig = sigfig)
  query_exp <- get_query_distribution(seu = seu, query = query, sigfig = sigfig)
  query_pct <- get_query_pct(seu = seu, query = query, by = "sample",
                             sigfig = sigfig)
  if("named_clusters" %in% colnames(seu@meta.data)) {
    get_by <- c("sample", "named_clusters")
  } else {
    get_by <- c("sample", "seurat_clusters")
  }
  query_cluster_pct <- get_query_pct(seu = seu, query = query, by = get_by,
                                     sigfig = sigfig)
  markers <- cbind(markers$gene, markers[, -grep("gene", colnames(markers))])
  colnames(markers)[1] <- "gene"
  message('Performing DEG analysis')
  DEG_list <- list()
  if(DEG_columns == "") {
    DEG_conditions <- int_columns
  }else{
    DEG_conditions <- unlist(strsplit(DEG_columns, split = ","))
  }
  for(condition in DEG_conditions) {
    message(paste0("Calculating DEGs for condition ", condition))
    condition_DEGs <- get_sc_markers(seu = seu, cond = condition, DEG = TRUE,
                                     verbose = verbose)
    DEG_list[[condition]] <- condition_DEGs
  }
  if(save_RDS){
  	message('Writing results to disk')
  	saveRDS(qc, file.path(output, paste0(project_name, ".qc.RDS")))
  	saveRDS(seu, file.path(output, paste0(project_name, ".seu.RDS")))
  	saveRDS(markers, file.path(output, paste0(project_name, ".markers.RDS")))
  	saveRDS(DEG_list, file.path(output, paste0(project_name, ".DEGs.RDS")))
  }
  final_results <- list()
  final_results$qc <- qc
  final_results$seu <- seu
  final_results$clusters_pct <- clusters_pct
  final_results$query_exp <- query_exp
  final_results$query_pct <- query_pct
  final_results$query_cluster_pct <- query_cluster_pct
  final_results$markers <- markers
  final_results$DEG_list <- DEG_list
  return(final_results)
}

#' write_seurat_report
#' Write integration HTML report
#' 
#' @param int integrated seurat object
#' @param output directory where report will be saved
#' @param name experiment name, will be used to build output file name
#' @param template_folder directory where template is located
#' @param source_folder htmlreportR source folder
#' @param int_columns factors present in experiment design
#'
#' @keywords preprocessing, write, report
#' 
#' @return nothing

write_seurat_report <- function(final_results, output = getwd(), name = NULL,
                                template_folder, source_folder = "none",
                                target_genes, int_columns, DEG_list = NULL,
                                cell_annotation = NULL){
  if(is.null(template_folder)) {
    stop("No template folder was provided.")
  }
  if(!file.exists(source_folder)) {
    stop(paste0("Source folder not found. Was ", source_folder))
  }
  if(any(is.null(final_results))) {
    stop("ERROR: final results object contains NULL fields. Analysis
       is not complete.")
  }
  template <- file.path(template_folder, "integration_template.txt")
  tmp_folder <- "tmp_lib"
  out_file <- file.path(output, "integration_report.html")
  container <- list(seu = final_results$seu, int_columns = int_columns,
                    DEG_list = final_results$DEG_list,
                    target_genes = target_genes,
                    clusters_pct = final_results$clusters_pct,
                    query_exp = final_results$query_exp,
                    query_pct = final_results$query_pct,
                    query_cluster_pct = query_cluster_pct,
                    markers = final_results$markers,
                    cell_annotation = cell_annotation)
  plotter <- htmlreportR:::htmlReport$new(title_doc = paste0("Single-Cell ", name, " report"), 
                            container = container, tmp_folder = tmp_folder,
                            src = source_folder)
  plotter$build(template)
  plotter$write_report(out_file)
  message(paste0("Report written in ", out_file))
}
