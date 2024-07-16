

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

#' do_qc
#' Perform Quality Control
#'
#' @param seu Seurat object
#' @param minqcfeats Min number of features for which a cell is selected
#' @param percentmt Max percentage of reads mapped to mitochondrial genes for which a cell is selected
#'
#' @keywords preprocessing, qc
#' 
#' @return Seurat object
do_qc <- function(seu, minqcfeats, percentmt){
  #### QC ####
  
  ##### Reads mapped to mitochondrial genes #####
  
  # In human ENSEMBL gene symbol annotations, mitochondrial genes are annotated starting with a MT- string
  
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  # We can define ribosomal proteins (their names start with RPS or RPL)
  
  seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  
  ##### Filtering out cells #####

  seu[['QC']] <- ifelse(seu@meta.data$Is_doublet == 'True','Doublet','Pass')
  seu[['QC']] <- ifelse(TRUE,'Pass','This should not happen') # provisional until I code the dublet detection stuff (see previous line)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_RNA < minqcfeats & seu@meta.data$QC == 'Pass','Low_nFeature',seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_RNA < minqcfeats & seu@meta.data$QC != 'Pass' & seu@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',seu@meta.data$QC,sep = ','),seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$percent.mt > percentmt & seu@meta.data$QC == 'Pass','High_MT',seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_RNA < minqcfeats & seu@meta.data$QC != 'Pass' & seu@meta.data$QC != 'High_MT',paste('High_MT',seu@meta.data$QC,sep = ','),seu@meta.data$QC)

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


##########################################################################


#' do_marker_gene_selection
#' Perform marker gene selection
#' TODO this function is harcoded - make proper variables
#'
#' @param seu Seurat object
#' @param out_path path where output will be saved. If not specified, it will
#' not be written to disk.
#' 
#' @keywords preprocessing, marker, gene
#' 
#' @return Nothing
do_marker_gene_selection <- function(seu, out_path = NULL){
  markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  if(!is.null(out_path)){ saveRDS(markers, paste0(out_path, ".markers.RDS"))}
  return(markers)
}

##########################################################################


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
    seu <- CreateSeuratObject(counts = d10x, project = sample, min.cells = 1,
                              min.features = 1)
    seu <- add_exp_design(seu = seu, name = sample, exp_design = exp_design)
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
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
#' @param anno_table Path to file containing cluster annotation, or loa
#'
#' @return Annotated seu object
#'

annotate_clusters <- function(seu, new_clusters = NULL ) {
  names(new_clusters) <- levels(seu)
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
  if(!is.null(seu@meta.data$named_clusters)) {
    cluster_idents <- "named_clusters"
  } else {
    warning("Seurat object clusters are not annotated.", immediate. = TRUE)
    cluster_idents <- "seurat_clusters"
  }
  if(DEG) {
    conds <- unique(seu[[cond]][[1]])
    if(length(conds) != 2) {
      stop(paste0("Error: Two groups must be supplied for DEG analysis.
                   Provided ", length(cond)))
    }
    return(list(idents = cond, conds = conds))
  } else {
    return(list(idents = cluster_idents, conds = NULL))
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
  }
  merged_df <- do.call(rbind, df_list)
  return(merged_df)
}

#' match_cell_types
#' `match_cell_types` takes a cluster-marker gene data frame and a cell type
#' marker file. It then looks for matches between the two and assigns a cell
#' type to each cluster of the data frame.
#'
#' @param markers_df Data frame of markers, clusters and p-values
#' @param anno_table Table of cell types and their associated markers
#' @param top Top markers by p-value to use in cell type assignment
#' @returns A markers data frame with a new column for cell type assigned to
#' cluster.

match_cell_types <- function(markers_df, anno_table, top = 20) {
  canon_types <- unique(anno_table$type)
  subset_list <- list()
  res <- list()
  pcols <- grep("p_val_adj", colnames(markers_df))
    if(length(pcols) > 1) {
      pvals <- unlist(markers_df[, pcols])
      pvals <- pvals[pvals != 0]
      min_pval <- min(pvals)
      markers_df$p_val_adj <- corto::fisherp(c(min(min_pval,
                                                   markers_df[[pcols[1]]]),
                                               min(min_pval,
                                                   markers_df[[pcols[2]]])))
    }
  fcols <- grep("log2FC", colnames(markers_df))
    if(length(fcols) > 1) {
      markers_df$avg_log2FC <- (markers_df[[pcols[1]]] + markers_df[[pcols[2]]]) / 2
    }
  for(cluster in unique(markers_df$cluster)) {
    subset <- markers_df[markers_df$cluster == cluster, ]
    subset <- subset[order(subset$p_val_adj), ]
    if(nrow(subset) > top) {
      subset <- subset[seq(1, top), ]
    }
    matches <- list()
    for(type in canon_types) {
      type_markers <- anno_table[anno_table$type == type, ]$marker
      # En vez de dividir, que sea una suma de coincidencias ponderada por media
      # de log2fc
      matches[[type]] <- sum(rownames(subset) %in% type_markers) /
                         length(type_markers)
    }
    if(max(unlist(matches)) == 0) {
      subset$cell_type = "Unknown"
    } else {
      cluster_match <- names(matches[which.max(unlist(matches))])
      subset$cell_type <- cluster_match
    }
    subset_list[[cluster]] <- subset
  }
  stats_table <- do.call(rbind, subset_list)
  stats_table <- stats_table[order(stats_table$cluster), ]
  stats_table$gene <- rownames(stats_table)
  columns <- colnames(stats_table)
  annotated_clusters <- stats_table[, columns %in% c("cluster", "cell_type"),
                                      drop = FALSE]
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
#' @param top An integer. Top conserved markers to select for each cluster.
#'            Default: 10
#' @returns A list containing one marker data frame per cluster.

get_sc_markers <- function(seu, cond = NULL, DEG = FALSE, top = 10) {
  cluster_idents <- get_marker_idents(seu = seu, cond = cond, DEG = DEG)
  Seurat::Idents(seu) <- cluster_idents$idents
  conds <- cluster_idents$conds
  clusters <- sort(unique(seu@meta.data[[cluster_idents$idents]]))
  clusters_markers <- list()
  for (i in seq(1, length(clusters))) {
    if (DEG) {
      message(paste0("Analysing factor ", i, "/", length(clusters)))
      # off-by-one correction because Seurat counts clusters from 0
      subset_seu <- subset_seurat(seu, "seurat_clusters", i - 1)
      markers <- Seurat::FindMarkers(subset_seu, ident.1 = conds[1],
                                     ident.2 = conds[2])
    } else {
      message(paste0("Analysing cluster ", i, "/", length(clusters)))
      markers <- Seurat::FindConservedMarkers(seu, ident.1 = clusters[i],
                                              grouping.var = cond,
                                              verbose = FALSE)[1:top, ]
    }
    clusters_markers[[as.character(clusters[i])]] <- markers
  }
  return(clusters_markers)
}

subset_seurat <- function(seu, column, value) {
  expr <- Seurat::FetchData(seu, vars = column)
  subset <- seu[, which(expr == value)]
  return(subset)
}

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

integrate_seurat <- function(seu, percentmt, ndims, resolution, hvgs, embeds,
                             minqcfeats, scalefactor, normalmethod,
                             clusters_annotation = NULL) {
  int <- do_qc(seu = seu, minqcfeats = minqcfeats, 
        percentmt = percentmt)
  int <- subset(int, subset = QC != 'High_MT,Low_nFeature')
  int <- Seurat::FindVariableFeatures(int, nfeatures = hvgs,
                                             selection.method = "vst")
  int <- Seurat::RunPCA(int)
  int <- harmony::RunHarmony(int, "sample",
                                    plot_convergence = FALSE)
  int <- Seurat::RunUMAP(int, dims = seq(1, ndims),
                        reduction = "harmony")
  int <- Seurat::FindNeighbors(int, dims = seq(1, ndims),
                        reduction = "harmony")
  int <- Seurat::FindClusters(int, resolution = 0.6)
  if(clusters_annotation != "") {
    message("Clusters annotation file provided. Annotating clusters")
    int <- annotate_clusters(seu = int, new_clusters = clusters_annotation[[2]])
    } else {
    message("Clusters annotation file not provided. Seurat object will
      be dynamically annotated")
    ## DEG, markers and de novo annotation will be included in function once
    ## they are stable
    }
  int <- SeuratObject::JoinLayers(int)
  return(int)
}

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
analyze_seurat <- function(raw_seu, out_path = NULL, minqcfeats, percentmt,
                           hvgs, ndims, resolution, dimreds_to_do,
                           embeds, integrate = FALSE){
  raw_seu <- do_qc(seu = raw_seu, minqcfeats = minqcfeats, 
                   percentmt = percentmt)
  seu <- subset(raw_seu, subset = QC != 'High_MT,Low_nFeature')
  seu <- raw_seu
  message('Finding variable features')
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = hvgs, verbose = FALSE)
  # dimreds_to_do: PCA/UMAP/tSNE
  message('Reducing dimensionality')
  seu <- do_dimred(seu = seu, ndims = ndims, dimreds = dimreds_to_do)   
  if(integrate) { # Harmony integration and remaining dimreds 
    message('Integration harmony and dimensionality reduction')
    seu <- harmony::RunHarmony(object = seu, group.by.vars = "sample",
                                verbose = FALSE)
    seu <- do_dimred(seu = seu, ndims = ndims, dimreds = c("tsne", "umap"),
                     reduction = embeds)
  }
  message('Clustering')
  seu <- do_clustering(seu = seu, ndims = ndims, resolution = resolution,
                       reduction = embeds)
  message('Selecting markers')
  markers <- SeuratObject::JoinLayers(seu)
  markers <- do_marker_gene_selection(seu = markers, out_path = out_path)
  if(!is.null(out_path)) saveRDS(seu, paste0(out_path, ".seu.RDS"))
  return(list(seu = seu, raw_seu = raw_seu, markers = markers))
}

##########################################################################


#' write_seurat_report
#' Write seurat HTML report
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
write_seurat_report <- function(all_seu = NULL, template, out_path,
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

#' write_integration_report
#' Write integration HTML report
#' 
#' @param int integrated seurat object
#' @param output_dir directory where report will be saved
#' @param name experiment name, will be used to build output file name
#' @param template_folder directory where template is located
#' @param source_folder htmlreportR source folder
#' @param int_columns factors present in experiment design
#'
#' @keywords preprocessing, write, report
#' 
#' @return nothing
write_integration_report <- function(int_seu, output_dir = getwd(),
                                     markers, template_folder, name = NULL,
                                     source_folder = "none", target_genes,
                                     int_columns, DEG_list = NULL,
                                     anno_table = NULL){
  if(is.null(template_folder)) {
    stop("No template folder was provided.")
  }
  if(!file.exists(source_folder)) {
    stop(paste0("Source folder not found. Was ", source_folder))
  }
  if(any(is.null(int_seu))) {
    stop("ERROR: comparison object contains NULL fields. Analysis
       is not complete.")
  }
  template <- file.path(template_folder, "integration_template.txt")
  tmp_folder <- "tmp_lib"
  out_file <- file.path(output_dir, "integration_report.html")
  container <- list(seu = int_seu, int_columns = int_columns,
                    DEG_list = DEG_list, target_genes = target_genes,
                    markers = markers, anno_table = anno_table)
  plotter <- htmlReport$new(title_doc = paste0("Single-Cell ", name, " report"), 
                            container = container, tmp_folder = tmp_folder,
                            src = source_folder)
  plotter$build(template)
  plotter$write_report(out_file)
}

#' format_markers
#' `format_markers` formats a marker-celltype table into a list
#' 
#' @param markers_list Table containing celltypes and their markers
#' 
#' @return A list with one element per cell type
read_and_format_markers <- function(path_to_markers) {
  markers_df <- read.table(path_to_markers, sep = "\t", header = FALSE,
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

