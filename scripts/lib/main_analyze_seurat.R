
#' main_analyze_seurat
#' `main_analyze_seurat` is the main seurat analysis function. Can be used
#' for integrative or non-integrative analysis.
#'
#' @param seu A seurat object.
#' @param minqcfeats An integer. Minimum features to consider a cell valid
#' @param percentmt A float. Maximun MT percentage to consider a cell valid
#' @param query A string vector. List of genes to explore in dataset
#' @param sigfig An integer. Significant figures to output
#' @param resolution An integer. Controls clustering granularity
#' @param p_adj_cutoff A float. Adjusted p-value cutoff by which to consider a maker
#' valid for cell type annotation
#' @param integrate A boolean.
#'   * `TRUE`: Integrate seurat object using the `harmony` package.
#'   * `FALSE` (the default): Do not perform integration.
#' @param cluster_annotation A data frame. Table to use to rename clusters.
#' @param cell_annotation A data frame. Table of markers to use for cell type
#' annotation
#' @param DEG_columns A string vector. Categories by which DEG analysis will be
#' performed
#' @param scalefactor An integer. Factor by which to scale data in normalisation
#' @param hvgs An integer. Number of highly-variable features to select
#' @param int_columns A string vector. Categories to consider in integrative
#' analysis.
#' @param normalmethod A string. Method to use in normalisation. Default is
#' "LogNormalize", the Seurat default.
#' @param ndims An integer. Target dimensions in dimensionality reduction.
#' @param verbose. A boolean.
#'   * `TRUE`: Prints progress bars and messages to closely monitor progress.
#'   * `FALSE` (the default): Print minimal progress messages.
#' @param output A string. Path to output directory.
#' @param save_RDS A boolean.
#'   * `TRUE`: Save output as an rds file.
#'   * `FALSE` (the default): Do not save output as an rds file.
#' @param reduce A boolean.
#'   * `TRUE`: Skip QC filtering. Intended for development and testing.
#'   * `FALSE` (the default): QC filtering will be performed.

main_analyze_seurat <- function(seu, minqcfeats, percentmt, query, sigfig = 2,
                           		  resolution, p_adj_cutoff = 5e-3,
                           		  integrate = FALSE, cluster_annotation = NULL,
                           		  cell_annotation = NULL, DEG_columns = NULL,
                           		  scalefactor = 10000, hvgs, int_columns = NULL,
                           		  normalmethod = "LogNormalize", ndims,
                                verbose = FALSE, output = getwd(),
                                save_RDS = FALSE, reduce = FALSE){
  qc <- tag_qc(seu = seu, minqcfeats = minqcfeats, percentmt = percentmt)
  colnames(qc@meta.data) <- tolower(colnames(qc@meta.data))
  if(!reduce) {
    seu <- subset(qc, subset = qc != 'High_MT,Low_nFeature')
  } else {
    message("Reduce argument is set to TRUE. Skipping QC subsetting")
    seu <- qc
  }
  message('Normalizing data')
  seu <- SeuratObject::JoinLayers(seu)
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
  seu <- Seurat::RunUMAP(object = seu, dims = seq(ndims),
                         reduction = reduction, verbose = verbose)
  seu <- Seurat::FindNeighbors(object = seu, dims = seq(1, ndims),
                               reduction = reduction, verbose = verbose)
  seu <- Seurat::FindClusters(seu, resolution = resolution, verbose = verbose)
  run_conserved <- ifelse(test = length(int_columns) == 1, no = FALSE,
                          yes = !has_exclusive_clusters(seu = seu,
                                                        cond = int_columns))
  if(run_conserved) {
    markers <- get_sc_markers(seu = seu, cond = int_columns, DEG = FALSE)
    markers <- collapse_markers(markers$markers)
  }else{
    markers <- Seurat::FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25,
                                      logfc.threshold = 0.25, verbose = verbose)
    rownames(markers) <- NULL
  }
  markers <- cbind(markers$gene, markers[, -grep("gene", colnames(markers))])
  colnames(markers)[1] <- "gene"
  if(!is.null(cluster_annotation)) {
  	message("Clusters annotation file provided. Annotating clusters.")
  	seu <- annotate_clusters(seu = seu, new_clusters = cluster_annotation$name)
  } else if(!is.null(cell_annotation)){
	  message(paste0("Clusters annotation file not provided. Dynamically ",
                   "annotating clusters."))
	  annotated_clusters <- match_cell_types(markers_df = markers,
                                           cell_annotation = cell_annotation,
                                           p_adj_cutoff = p_adj_cutoff)
	  markers <- annotated_clusters$summary
	  seu <- annotate_clusters(seu, annotated_clusters$cell_types)
  } else {
  	warning("No data provided for cluster annotation.")
  }
  message("Extracting expression quality metrics.")
  sample_qc_pct <- get_qc_pct(seu, by = "sample")
  message("Extracting query expression metrics. This might take a while.")
  clusters_pct <- get_clusters_distribution(seu = seu, sigfig = sigfig)
  if(!is.null(query)) {
    query_data <- analyze_query(seu = seu, query = query, sigfig = sigfig)
  } else {
    query_data <- NULL
  }
  if(!is.null(DEG_columns)) {
    message('Performing DEG analysis.')
    DEG_conditions <- unlist(strsplit(DEG_columns, split = ","))
    if(length(DEG_conditions) == 0) {
      DEG_conditions <- int_columns
    }
    DEG_list <- vector(mode = "list", length = length(DEG_conditions))
    names(DEG_list) <- DEG_conditions
    DEG_meta <- DEG_conditions
    for(condition in DEG_conditions) {
      message(paste0("Calculating DEGs for condition ", condition, "."))
      condition_DEGs <- get_sc_markers(seu = seu, cond = condition, DEG = TRUE,
                                       verbose = verbose)
      DEG_list[[condition]] <- condition_DEGs
    }
    subset_DEGs <- NULL
    subset_seu <- NULL
    if(length(int_columns) == 2) {
      message(paste0("Analysing DEGs by subgroups. Subsetting by condition: ",
              DEG_conditions[1], ", analyzing effects of ", DEG_conditions[2]))
      subset_DEGs <- vector(mode = "list", length = 2)
      condition_values <- unique(seu@meta.data[[DEG_conditions[1]]])
      names(subset_DEGs) <- condition_values
      subset_seu <- subset_DEGs
      for(value in condition_values) {
        subset_seu[[value]] <- subset_seurat(seu, DEG_conditions[1], value)
        subset_DEGs[[value]] <- get_sc_markers(seu = subset_seu[[value]],
                                  cond = DEG_conditions[2], DEG = TRUE,
                                  verbose = verbose)
      }
    }
  } else {
    DEG_list <- NULL
  }
  final_results <- list()
  final_results$qc <- qc
  final_results$seu <- seu
  final_results$sample_qc_pct <- sample_qc_pct
  final_results$clusters_pct <- clusters_pct
  final_results$query_exp <- query_data$query_exp
  final_results$query_pct <- query_data$query_pct
  final_results$query_cluster_pct <- query_data$query_cluster_pct
  final_results$markers <- markers
  final_results$DEG_list <- DEG_list
  final_results$subset_seu <- subset_seu
  final_results$subset_DEGs <- subset_DEGs
  if(save_RDS){
    message('Writing results to disk.')
    saveRDS(final_results, file.path(output, paste0(seu@project.name, ".final_results.rds")))
  }
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
                                target_genes = NULL, int_columns = NULL,
                                DEG_list = NULL, cell_annotation = NULL,
                                template = NULL, out_name = NULL){
  if(is.null(template_folder)) {
    stop("No template folder was provided.")
  }
  if(!file.exists(source_folder)) {
    stop(paste0("Source folder not found. Was ", source_folder, " ."))
  }
  if(any(is.null(final_results))) {
    stop("ERROR: final results object contains NULL fields. Analysis
       is not complete.")
  }
  if(is.null(template)) {
    stop("Please specify a template to render")
  }
  template <- file.path(template_folder, template)
  if(!file.exists(template)) {
    stop("Specified template does not exist in template folder.")
  }
  tmp_folder <- "tmp_lib"
  out_file <- file.path(output, paste0(name, "_", out_name))
  container <- list(seu = final_results$seu, int_columns = int_columns,
                    DEG_list = final_results$DEG_list,
                    marker_meta = final_results$marker_meta,
                    subset_seu = final_results$subset_seu,
                    subset_DEGs = final_results$subset_DEGs,
                    target_genes = target_genes,
                    sample_qc_pct = final_results$sample_qc_pct,
                    clusters_pct = final_results$clusters_pct,
                    query_exp = final_results$query_exp,
                    query_pct = final_results$query_pct,
                    query_cluster_pct = final_results$query_cluster_pct,
                    markers = final_results$markers,
                    cell_annotation = cell_annotation)
  plotter <- htmlreportR:::htmlReport$new(title_doc = paste0("Single-Cell ", name, " report"), 
                            container = container, tmp_folder = tmp_folder,
                            src = source_folder)
  plotter$build(template)
  plotter$write_report(out_file)
  message(paste0("Report written in ", out_file))
}
