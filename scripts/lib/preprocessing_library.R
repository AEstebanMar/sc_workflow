

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
  seu[['qc']] <- ifelse(seu@meta.data$nfeature_rna < minqcfeats,
                        'Low_nFeature', 'Pass')
  seu[['qc']] <- ifelse(seu@meta.data$nfeature_rna < minqcfeats &
                        seu@meta.data$qc != 'Pass' &
                        seu@meta.data$qc != 'Low_nFeature',
                        paste('Low_nFeature', seu@meta.data$qc, sep = ','),
                              seu@meta.data$qc)
  seu[['qc']] <- ifelse(seu@meta.data$percent.mt > percentmt &
                        seu@meta.data$qc == 'Pass','High_MT', seu@meta.data$qc)
  seu[['qc']] <- ifelse(seu@meta.data$nfeature_rna < minqcfeats &
                        seu@meta.data$qc != 'Pass' &
                        seu@meta.data$qc != 'High_MT',
                        paste('High_MT', seu@meta.data$qc ,sep = ','),
                              seu@meta.data$qc)
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
  seu <- Seurat::RenameIdents(seu, new_clusters)
  seu@meta.data$named_clusters <- Seurat::Idents(seu)
  return(seu)
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
  merged_df <- cbind(merged_df$gene, merged_df[, colnames(merged_df) != "gene"])
  colnames(merged_df)[1] <- "gene"
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
    pvals <- pvals[pvals != 0]
    min_pval <- min(pvals) # Small correction for machine-zeroes
    int_p_val <- rep(0, nrow(markers_df))
    p_vals_1 <- markers_df[[pcols[1]]]
    p_vals_2 <- markers_df[[pcols[2]]]
    for(i in seq(1, length(int_p_val))) {
      int_p_val[i] <- corto::fisherp(c(max(min_pval, p_vals_1[i]),
                                       max(min_pval, p_vals_2[i])))
    }
    markers_df$p_val_adj <- int_p_val
  } else {
    colnames(markers_df)[pcols] <- "p_val_adj"
  }
  markers_df <- markers_df[markers_df$p_val_adj <= p_adj_cutoff, ]
  fcols <- grep("log2FC", colnames(markers_df))
  if(any(is.na(markers_df[, fcols]))) {
    warning("WARNING: NAs detected in marker log2FC. Coercing to 0.",
            immediate. = TRUE)
    markers_df[, fcols][is.na(markers_df[, fcols])] <- 0
  }
  if(length(fcols) > 1) {
    markers_df$avg_log2FC <- (markers_df[[fcols[1]]] +
                              markers_df[[fcols[2]]]) / 2
  } else {
    colnames(markers_df)[fcols] <- "avg_log2FC"
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
      cluster_match <- "Unknown"
    } else {
      cluster_match <- names(scores[which.max(unlist(scores))])
    }
    subset$cell_type <- paste0(subset$cluster, ". ", cluster_match)
    subset_list[[as.numeric(cluster) + 1]] <- subset
  }
  stats_table <- do.call(rbind, subset_list)
  stats_table <- stats_table[order(stats_table$cluster), ]
  columns <- colnames(stats_table)
  anno_types <- strsplit(stats_table$cell_type, "\\. ")
  anno_types <- sapply(anno_types, `[`, 2)
  types <- unique(anno_types)
  for(type in types) {
    matches <- which(anno_types == type)
    type_clusters <- stats_table$cluster[matches]
    if(length(unique(type_clusters)) > 1) {
      index <- as.integer(as.factor(type_clusters))
      dupe_cluster <- paste0(stats_table$cell_type[matches],
                             " (", letters[index], ")")
      stats_table[matches, ]$cell_type <- dupe_cluster
    }
  }
  sum_columns <- c("gene", "p_val_adj", "avg_log2FC", "cluster", "cell_type")
  res <- list(stats_table = stats_table,
              cell_types = unique(stats_table$cell_type),
              summary = stats_table[, sum_columns])
  return(res)
}

#' get_sc_markers
#' `get_sc_markers` performs differential expression analysis on OR selects
#' conserver markers from a seurat object.
#'
#' @param seu Seurat object to analyse.
#' @param DEG A boolean.
#'   * `TRUE`: Function will calculate differentally expressed genes.
#'   * `FALSE` (the default): Function will calculate cluster markers.
#' @param cond A string. Condition by which to perform DEG analysis, or by which
#' to group data to find conserved markers.
#' @param verbose A boolean. Will be passed to Seurat function calls.
#' @returns A list containing one marker or DEG data frame per cluster, plus
#' an additional one for global DEGs if performing differential analysis.

get_sc_markers <- function(seu, cond = NULL, DEG = FALSE, verbose = FALSE) {
  Seurat::Idents(seu) <- "seurat_clusters"
  conds <- unique(seu@meta.data[[cond]])
  clusters <- sort(unique(Seurat::Idents(seu)))
  cluster_markers <- list()
  for (i in seq(1, length(clusters))) {
    message(paste0("Analysing cluster ", i, "/", length(clusters)))
    if(DEG) {
      # off-by-one correction because Seurat counts clusters from 0
      subset_seu <- subset_seurat(seu, "seurat_clusters", i - 1)
      meta <- subset_seu@meta.data[[cond]]
      ncells <- c(sum(meta==conds[1]), sum(meta==conds[2]))
      if(any(ncells < 3)) {
        warning(paste0('Cluster ', i, ' contains less than three cells for',
                        ' condition \'', conds[which(ncells < 3)],
                        '\'. Skipping DEG analysis', collapse = ""),
                immediate. = TRUE)
        markers <- data.frame(FALSE)
      } else {
        Seurat::Idents(subset_seu) <- cond  
        markers <- Seurat::FindMarkers(subset_seu, ident.1 = conds[1],
                                       ident.2 = conds[2], verbose = verbose)
        markers$gene <- rownames(markers)
      }
      } else {
      markers <- Seurat::FindConservedMarkers(seu, ident.1 = clusters[i],
                                              grouping.var = cond,
                                              verbose = verbose)
    }
    nums <- sapply(markers, is.numeric)
    markers[nums] <- lapply(markers[nums], signif, 2)
    cluster_markers[[as.character(clusters[i])]] <- markers
  }
  if(DEG) {
    message("Calculating global DEGs")
    Seurat::Idents(seu) <- cond
    global_markers <- Seurat::FindMarkers(seu, ident.1 = conds[1],
                                ident.2 = conds[2], verbose = verbose)
    global_markers$gene <- rownames(global_markers)
    nums <- sapply(global_markers, is.numeric)
    global_markers[nums] <- lapply(global_markers[nums], signif, 2)
    cluster_markers[["global"]] <- global_markers
  }
  if(!is.null(seu@meta.data$named_clusters)) {
    metadata <- unique(seu@meta.data[, c("seurat_clusters", "named_clusters")])
    named_clusters <- metadata[order(metadata$seurat_clusters), ]$named_clusters
    names(cluster_markers) <- c(as.character(named_clusters), "global")
  }
  return(cluster_markers)
}

#' get_clusters_distribution
#' `get_clusters_distribution` calculates the percentage of cells that make up
#' each cluster for each different sample in a seurat object. If clusters are
#' annotated, it will show cell types instead of cluster number.
#'
#' @param seu Clustered seurat object.
#' @param sigfig Significant figure cutoff
#' @returns A data frame with cell type distribution in each sample.

get_clusters_distribution <- function(seu, sigfig = 3) {
  clusters_column <- ifelse("named_clusters" %in% colnames(seu@meta.data), 
                            "named_clusters", "seurat_clusters")
  clusters_table <- table(seu@meta.data[, c("sample", clusters_column)])
  percent_table <- signif(clusters_table/rowSums(clusters_table)*100, sigfig)
  percent_table <- as.data.frame.matrix(percent_table)
  return(percent_table)
}

#' get_query_distribution
#' `get_query_distribution` builds a table of query genes expression levels
#' across all samples of a Seurat object
#'
#' @param seu Seurat object
#' @param query Vector of query genes whose expression to analyse.
#' @param sigfig Significant figure cutoff
#' @return A data frame with expression levels for query genes in each sample.

get_query_distribution <- function(seu, query, sigfig = 3) {
  genes <- SeuratObject::FetchData(seu, query)
  genes <- cbind(seu@meta.data$sample, genes)
  colnames(genes)[1] <- "sample"
  gene_distribution <- aggregate(genes[, -1], list(genes$sample), FUN = sum)
  gene_distribution[, -1] <- signif(gene_distribution[, -1], sigfig)
  rownames(gene_distribution) <- gene_distribution[, 1]
  gene_distribution <- gene_distribution[, -1, drop = FALSE]
  # In the case where query vector is of length one, its names are dropped.
  # Therefore, we need to set them forcefully.
  colnames(gene_distribution) <- query
  return(gene_distribution)
}

#' get_query_pct
#' `get_query_pct` gets the percentage of cells in each sample of a seurat
#' object which expresses genes specified in a list of queries.
#'
#' @param seu Seurat object
#' @param query Vector of query genes whose expression to analyse.
#' @param sigfig Significant figure cutoff
#' @param assay Seurat assay from which to extract data. Default is "RNA",
#' the default assay.
#' @param layer Layer of Seurat object from which to extract data. Default is
#' "counts", normalised assay data.
#' @return A data frame with expression levels for query genes in each sample.

get_query_pct <- function(seu, query, sigfig = 2, assay = "RNA",
                          layer = "counts") {
  pct_list <- list()
  for(sample in unique(seu@meta.data$sample)) {
    subset_seu <- subset_seurat(seu, "sample", sample)
    genes <- SeuratObject::GetAssayData(seu, assay = assay,
                                                    layer = layer)
    missing <- !(query %in% rownames(genes))
    if(any(missing)) {
      warning(paste0("Query genes ", paste0(query[missing],
                                            collapse = ", "),
                     " not present in seurat object."), immediate. = TRUE)
      query <- query[!missing]
    }
    queries <- genes[query, ]
    if(is.vector(queries)) {
      pct <- mean(queries)
      names(pct) <- query
    } else {
      pct <- rowMeans(queries)
    }
    pct_list[[sample]] <- pct
  }
  res <- do.call(rbind, pct_list) * 100
  res <- signif(res, sigfig)
  return(res)
}

#' has_exclusive_clusters
#' `has_exclusive_clusters` checks if seurat object contains clusters with
#' less than three members (cells) for any category in provided experimental
#' condition.
#'
#' @param seu Seurat object
#' @param cond Experimental condition
#' @returns A boolean. TRUE if exclusive clusters are found, FALSE otherwise.

has_exclusive_clusters <- function(seu, cond) {
  meta <- seu@meta.data[, c(cond, "seurat_clusters")]
  groups <- unique(meta[[cond]])
  clusters <- unique(meta[["seurat_clusters"]])
  pairs <- expand.grid(groups, clusters)
  sum_matches <- vector(mode = "integer")
  for(pair in seq(1, nrow(pairs))) {
    matches <- apply(meta, 1, function(x) x == pairs[pair, ])
    sum_matches[pair] <- sum(colSums(matches) == 2)
  }
  if(any(sum_matches < 3)) {
    mismatch <- pairs[which(sum_matches < 3), ]
    warning('One or more clusters contain less than three cells for one or ',
            'more categories. Affected pair(s): ',
            paste(apply(mismatch, 1, paste, collapse = "-"), collapse = ", "),
            ". \nDefaulting to general marker analysis.")
    res <- TRUE
  }else{
    res <- FALSE
  }
  return(res)
}

#' subset_seurat
#' `subset_seurat` subsets a seurat object. Documentation in progress.

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

