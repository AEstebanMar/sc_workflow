

#' read_sc_counts
#' Create seurat object from cellranger counts
#'
#' @param name sample name
#' @param input path to cellranger counts
#' @param mincells min number of cells for which a feature is recorded
#' @param minfeats min number of features for which a cell is recorded
#' 
#' @return Seurat object
read_input <- function(name, input, mincells, minfeats){
  mtx <- Read10X(input)
  seu <- CreateSeuratObject(counts = mtx, project = name, min.cells = mincells, 
                            min.features = minfeats)
  return(seu)
}

#' tag_gc
#' Perform Quality Control
#'
#' @param seu Seurat object
#' @param minqcfeats Min number of features for which a cell is selected
#' @param percentmt Max percentage of reads mapped to mitochondrial genes for which a cell is selected
#'
#' @keywords preprocessing, qc
#' 
#' @return Seurat object
tag_qc <- function(seu, minqcfeats, percentmt){
  
  seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["percent.rb"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]")
  # TO DO: doublet detection
  # seu[['QC']] <- ifelse(seu@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
  # seu[['QC']] <- ifelse(TRUE, 'Pass', 'This should not happen')
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_RNA < minqcfeats &
                        seu@meta.data$QC == 'Pass','Low_nFeature',
                        seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_RNA < minqcfeats &
                        seu@meta.data$QC != 'Pass' &
                        seu@meta.data$QC != 'Low_nFeature',
                        paste('Low_nFeature', seu@meta.data$QC, sep = ','),
                              seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$percent.mt > percentmt &
                        seu@meta.data$QC == 'Pass','High_MT',seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_RNA < minqcfeats &
                        seu@meta.data$QC != 'Pass' &
                        seu@meta.data$QC != 'High_MT',
                        paste('High_MT', seu@meta.data$QC,sep = ','),
                              seu@meta.data$QC)
  return(seu)
}


##########################################################################


#' do_dimred
#' Perform linear (PCA) and non-linear (UMAP/tSNE) dimensionality reduction
#'
#' @param seu Seurat object
#' @param ndims Number of PC to be used for UMAP / tSNE
#' @param dimreds character vector with the dimensional reductions to perform. E.g. c("pca", "tsne", "umap")
#' @param reduction Dimensional reduction to use for UMAP /tSNE. "pca" if no integration, or "harmony" if integration
#'
#' @keywords preprocessing, dimensionality, reduction, PCA, UMAP, tSNE
#' 
#' @return Seurat object
do_dimred <- function(seu, ndims, dimreds, reduction = "pca"){
  if ("pca" %in% dimreds){
    seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  }
  if ("umap" %in% dimreds){
    seu <- RunUMAP(seu, dims = 1:ndims, reduction = reduction)
  }
  if ("tsne" %in% dimreds){
    seu <- RunTSNE(seu, dims = 1:ndims, reduction = reduction)
  }
  return(seu)
}


##########################################################################


#' do_clustering
#' Perform clustering of cells
#'
#' @param seu Seurat object
#' @param ndims Number of PC to be used for clustering
#' @param resolution Granularity of the downstream clustering (higher values -> greater number of clusters)
#' @param reduction Dimensional reduction to use for clustering. "pca" if no integration, or "harmony" if integration
#' 
#' @keywords preprocessing, clustering
#' 
#' @return Seurat object
do_clustering <- function(seu, ndims, resolution, reduction){
  seu <- FindNeighbors(seu, dims = 1:ndims, reduction = reduction)
  seu <- FindClusters(seu, resolution = resolution)
  return(seu)
}

#' add_exp_design
#' Add experimental condition to Seurat metadata
#'
#' @param seu Seurat object
#' @param name Sample name
#' @param exp_design Experiment design table in CSV format
#' 
#' @keywords preprocessing, subsetting, integration
#' 
#' @return Seurat object with the experimental conditions added as metadata
add_exp_design <- function(seu, name, exp_design){
  exp_design <- as.list(exp_design[exp_design$sample == name,])
  for (i in names(exp_design)){
    seu@meta.data[[i]] <- c(rep(exp_design[[i]], nrow(seu@meta.data)))
  }
  return(seu)
}


##########################################################################

#' merge_seurat
#' `merge_seurat` loads single-cell count matrices and creates a merged
#' seurat object.
#'
#' @param exp_cond Experimental condition
#' @param samples Vector of samples with that experimental condition
#' @param exp_design Experiment design table in CSV format
#' @param count_path Directory with count results
#' 
#' @keywords preprocessing, merging, integration
#' 
#' @return Merged Seurat object
merge_seurat <- function(project_name, samples, exp_design, count_path,
                         suffix=''){
  full_paths <- Sys.glob(paste(count_path, suffix, sep = "/"))
  seu.list <- sapply(samples, function(sample){
    sample_path <- grep(sample, full_paths, value = TRUE)
    d10x <- Seurat::Read10X(sample_path)
    seu <- Seurat::CreateSeuratObject(counts = d10x, project = sample, min.cells = 1,
                              min.features = 1)
    seu <- add_exp_design(seu = seu, name = sample, exp_design = exp_design)
    seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
    seu <- subset(seu, subset = nFeature_RNA > 500 & nCount_RNA > 1000
                                & percent.mt < 20)
    return(seu)
    })
  merged_seu <- merge(seu.list[[1]], y = seu.list[-1], add.cell.ids = samples,
                      project = project_name)
  return(merged_seu)
}

#' annotate_clusters
#' `annotate_clusters` renames seurat clusters according to dictionary
#'
#' @param seu Non-annotated seu with markers
#' @param new_clusters Vector of names to assign to clusters.
#'
#' @return Annotated seu object
#'

annotate_clusters <- function(seu, new_clusters = NULL ) {
  names(new_clusters) <- levels(seu)
  types <- unique(new_clusters)
  for(type in types) {
    matches <- new_clusters %in% type
    marked_dupes <- paste0(new_clusters, " (", letters[cumsum(matches)], ")")
    new_clusters[matches] <- marked_dupes[matches]
  }
  seu <- Seurat::RenameIdents(seu, new_clusters)
  seu@meta.data$named_clusters <- Seurat::Idents(seu)
  return(seu)
}

#' get_marker_idents
#' `get_marker_idents` performs some checks in input seurat object metadata. It
#' checks if it is clustered, if clusters are annotated, and, if DEG mode is
#' active, it checks whether provided condition allows DEG analysis.
#'
#' @param seu Seurat object on which DEG analysis will be performed.
#' @param DEG A boolean.
#'   * `TRUE`: Function will calculate differentally expressed genes.
#'   * `FALSE` (the default): Function will calculate conserved markers.
#' @param cond A string. Condition by which to perform DEG analysis.
#' @returns A list containing items `idents`, with new seurat identifiers,
#' and item `conds`, which, if DEG mode is active, contains unique values
#' of specified column of seurat metadata. If not active, it is set to NULL.
#' '

get_marker_idents <- function(seu, cond, DEG) {
  if(is.null(seu@meta.data$seurat_clusters)) {
    stop("Error: Seurat object contains no clusters. Analysis impossible.")
  }
  if(!cond %in% names(seu@meta.data)) {
    stop("Specified condition does not exist in seurat object metadata.
          DEG analysis impossible.")
  }
  if(DEG) {
    conds <- unique(seu[[cond]][[1]])
    if(length(conds) != 2) {
      stop(paste0("Error: Two groups must be supplied for DEG analysis.
                   Provided ", length(cond)))
    }
    return(list(idents = cond, conds = conds))
  } else {
    return(list(idents = "seurat_clusters", conds = NULL))
  }
}

#' collapse_markers
#' `collapse_markers` takes list of marker gene data frames and collapses it
#' into a cluster-markers data frame.
#'
#' @param marker_list A list containing marker gene data frames.
#' @returns A data frame. Column `cluster` contains element names of original
#' list (seurat clusters) and column `genes` contains, for each row, the top
#' markers of that element from the original list separated by commas.

collapse_markers <- function(markers_list) {
  df_list <- list()
  for(i in seq(1, length(markers_list))) {
    df_list[[i]] <- as.data.frame(markers_list[[i]])
    df_list[[i]]$cluster <- i
    df_list[[i]]$gene <- rownames(df_list[[i]])
    rownames(df_list[[i]]) <- NULL
  }
  merged_df <- do.call(plyr::rbind.fill, df_list)
  merged_df
  return(merged_df)
}

#' match_cell_types
#' `match_cell_types` takes a cluster-marker gene data frame and a cell type
#' marker file. It then looks for matches between the two and assigns a cell
#' type to each cluster of the data frame.
#'
#' @param markers_df Data frame of markers, clusters and p-values
#' @param cell_annotation Table of cell types and their associated markers
#' @param top Top markers by p-value to use in cell type assignment
#' @returns A markers data frame with a new column for cell type assigned to
#' cluster.

match_cell_types <- function(markers_df, cell_annotation, p_adj_cutoff = 1e-5) {
  canon_types <- unique(cell_annotation$type)
  subset_list <- list()
  res <- list()
  pcols <- grep("p_val_adj", colnames(markers_df))
  if(any(is.na(markers_df[, pcols]))) {
    warning("WARNING: NAs detected in marker p-values. Coercing to 1.",
            immediate. = TRUE)
    markers_df[, pcols][is.na(markers_df[, pcols])] <- 1
  }
  if(length(pcols) > 1) {
    pvals <- unlist(markers_df[, pcols])
    pvals <- pvals[pvals != 0 & !is.na(pvals)]
    min_pval <- min(pvals)
    markers_df$p_val_adj <- corto::fisherp(c(min(min_pval,
                                                 markers_df[[pcols[1]]]),
                                             min(min_pval,
                                                 markers_df[[pcols[2]]])))
  }
  markers_df <- markers_df[markers_df$p_val_adj < p_adj_cutoff, ]
  fcols <- grep("log2FC", colnames(markers_df))
  if(length(fcols) > 1) {
    markers_df$avg_log2FC <- (markers_df[[fcols[1]]] +
                              markers_df[[fcols[2]]]) / 2
  }
  for(cluster in unique(markers_df$cluster)) {
    subset <- markers_df[markers_df$cluster == cluster, ]
    subset <- subset[order(subset$p_val_adj), ]
    max_log2FC <- max(subset$avg_log2FC)
    scores <- list()
    for(type in canon_types) {
      type_markers <- cell_annotation[cell_annotation$type == type, ]$marker
      found_markers <- which(subset$gene %in% type_markers)
      scores[[type]] <- sum(subset$avg_log2FC[found_markers] / max_log2FC)
    }
    if(max(unlist(scores)) == 0 || is.na(max(unlist(scores)))) {
      subset$cell_type = "Unknown"
    } else {
      cluster_match <- names(scores[which.max(unlist(scores))])
      subset$cell_type <- cluster_match
    }
    subset_list[[cluster]] <- subset
  }
  stats_table <- do.call(rbind, subset_list)
  stats_table <- stats_table[order(stats_table$cluster), ]
  stats_table$gene <- rownames(stats_table)
  columns <- colnames(stats_table)
  annotated_clusters <- stats_table[, c("cluster", "cell_type")]
  rownames(annotated_clusters) <- NULL
  res <- list(stats_table = stats_table,
              cell_types = unique(annotated_clusters)$cell_type)
  return(res)
}

#' get_sc_markers
#' `get_sc_markers` performs differential expression analysis on OR selects
#' conserver markers from a seurat object.
#'
#' @inheritParams get_marker_idents
#' @returns A list containing one marker data frame per cluster.

get_sc_markers <- function(seu, cond = NULL, DEG = FALSE) {
  cluster_idents <- get_marker_idents(seu = seu, cond = cond, DEG = DEG)
  Seurat::Idents(seu) <- cluster_idents$idents
  conds <- cluster_idents$conds
  clusters <- sort(unique(seu@meta.data[[cluster_idents$idents]]))
  cluster_markers <- list()
  for (i in seq(1, length(clusters))) {
    if (DEG) {
      message(paste0("Analysing factor ", i, "/", length(clusters)))
      # off-by-one correction because Seurat counts clusters from 0
      subset_seu <- subset_seurat(seu, "seurat_clusters", i - 1)
      markers <- Seurat::FindMarkers(subset_seu, ident.1 = conds[1],
                                     ident.2 = conds[2], verbose = FALSE)
    } else {
      message(paste0("Analysing cluster ", i, "/", length(clusters)))
      markers <- Seurat::FindConservedMarkers(seu, ident.1 = clusters[i],
                                              grouping.var = cond,
                                              verbose = FALSE)
      markers$gene <- rownames(markers)
      rownames(markers) <- NULL
      markers <- cbind(markers$gene, markers[,
                                              -grep("gene", colnames(markers))])
    }
    cluster_markers[[as.character(clusters[i])]] <- markers
  }
  return(cluster_markers)
}

subset_seurat <- function(seu, column, value) {
  expr <- Seurat::FetchData(seu, vars = column)
  subset <- seu[, which(expr == value)]
  return(subset)
}

#' write_sergio_report
#' Write sergio's seurat HTML report
#' 
#' @param name sample name
#' @param expermient experiment name
#' @param template Rmd template
#' @param outdir output directory
#' @param intermediate_files directory for saving intermediate files in case pandoc fails
#' @param minqcfeats min number of features for which a cell is selected
#' @param percentmt max percentage of reads mapped to mitochondrial genes for which a cell is selected
#' @param hvgs Number of HVGs to be selected
#' @param resolution Granularity of the downstream clustering (higher values -> greater number of clusters)
#' @param all_seu NULL if creating an individual report (daemon 3a). A list of 2 Seurat objects and a matrix of markers if creating a general report (daemon 3b)
#' 
#' @keywords preprocessing, write, report
#' 
#' @return nothing

write_sergio_report <- function(all_seu = NULL, template, out_path,
                                intermediate_files, minqcfeats,
                                percentmt, hvgs, resolution){
  int_files <- file.path(out_path, intermediate_files)
  if (!file.exists(int_files)) dir.create(int_files)
  if (is.null(all_seu)){
    seu <- readRDS(paste0(out_path, ".seu.RDS"))
    before.seu <- readRDS(paste0(out_path, ".before.seu.RDS"))
    markers <- readRDS(paste0(out_path, ".markers.RDS"))
  } else {
    seu <- all_seu[[1]]
    before.seu <- all_seu[[2]]
    markers <- all_seu[[3]]
  }  
  rmarkdown::render(template, clean = TRUE, intermediates_dir = int_files,
                    output_file = paste0(out_path,
                                      "_preprocessing_report.html"))
}

#' read_and_format_targets
#' `read_and_format_targets` formats a marker-celltype table into a list
#' 
#' @param file Path to target gene file
#' 
#' @return A list with one element per cell type

read_and_format_targets <- function(file) {
  markers_df <- read.table(file, sep = "\t", header = FALSE,
                           stringsAsFactors = FALSE)
  cell_types <- markers_df[, 1]
  markers <- strsplit(markers_df[, 2], ",")
  names(markers) <- cell_types
  return(markers)
}


##########################################################################


#' extract_metadata
#' Extract metadata dataframe from Seurat objects
#' 
#' @param seu Seurat object / list of Seurat objects
#' 
#' @keywords preprocessing, report, metadata
#' 
#' @return Dataframe with metadata
extract_metadata <- function(seu){
  if (!is.list(seu)){
    seu <- seu[[]]
  } else {
    seu <- lapply(seu, "[[")
    seu <- do.call(rbind, seu)
    }
return(seu)
}

##########################################################################

#' make_vln
#' Make Violin plot
#' 
#' @param seu Seurat object / list of Seurat objects
#' @param feature metadata feature to plot
#' 
#' @keywords preprocessing, report, plot, violin
#' 
#' @return nothing
make_vln <- function(seu, feature){
  seu <- extract_metadata(seu)
  seu <- seu[, c("orig.ident", feature)]
  colnames(seu)[2] <- "values"
  ggplot() + 
    geom_point(seu,
               mapping = aes(orig.ident, values, fill = orig.ident),
               shape = 21,
               colour = "white",
               size = 2,
               stroke = 0.5,
               position = "jitter",
               alpha = 0.3) +
    geom_violin(seu,
                mapping = aes(orig.ident, values, fill = orig.ident),
                width = 0.5,
                color = "black") +
    xlab(NULL) +
    ylab(feature) +
    labs(fill = NULL) +
    theme_bw()
}

##########################################################################


#' ensure_list
#' Makes sure you have list of Seurat objects (even if you have only one)
#' 
#' @param seu Seurat object / list of Seurat objects
#' 
#' @keywords preprocessing, report, list
#' 
#' @return List of Seurat objects
ensure_list <- function(seu){
  if (!is.list(seu)){
    seu <- list(seu)
  }
  return(seu)
}

##########################################################################

