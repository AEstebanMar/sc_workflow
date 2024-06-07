# Sergio Alías, 20230606
# Last modified 20231228

##########################################################################
########################## PRE-PROCESSING LIBRARY ########################
##########################################################################


#' append_message
#' Append experiments names with a Hello World message. Useful for testing :D
#'
#' @param name expermient name
#' @param message message to append. Default: "Hello World"
#' 
#' @keywords test
#' 
#' @return character vector
append_message <- function(name, message = "Hello World") {paste(name, message)}


##########################################################################


#' read_input
#' Read Cellranger input into Seurat
#'
#' @param name sample name
#' @param input directory with the single-cell data
#' @param mincells min number of cells for which a feature is recorded
#' @param minfeats min number of features for which a cell is recorded
#' 
#' @keywords preprocessing, input
#' 
#' @return Seurat object
read_input <- function(name, input, mincells, minfeats){
  mtx <- Read10X(input)
  seu <- CreateSeuratObject(counts = mtx,
                            project = name, 
                            assay = "scRNAseq", 
                            min.cells = mincells, 
                            min.features = minfeats
                            )
  mtx <- NULL
  return(seu)
}


##########################################################################


#' do_qc
#' Perform Quality Control
#'
#' @param name sample name
#' @param expermient experiment name
#' @param seu Seurat object
#' @param minqcfeats Min number of features for which a cell is selected
#' @param percentmt Max percentage of reads mapped to mitochondrial genes for which a cell is selected
#' @param save_before_seu Whether to save Seurat object before QC
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

  # seu[['QC']] <- ifelse(seu@meta.data$Is_doublet == 'True','Doublet','Pass')
  seu[['QC']] <- ifelse(TRUE,'Pass','This should not happen') # provisional until I code the dublet detection stuff (see previous line)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_scRNAseq < minqcfeats & seu@meta.data$QC == 'Pass','Low_nFeature',seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_scRNAseq < minqcfeats & seu@meta.data$QC != 'Pass' & seu@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',seu@meta.data$QC,sep = ','),seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$percent.mt > percentmt & seu@meta.data$QC == 'Pass','High_MT',seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_scRNAseq < minqcfeats & seu@meta.data$QC != 'Pass' & seu@meta.data$QC != 'High_MT',paste('High_MT',seu@meta.data$QC,sep = ','),seu@meta.data$QC)

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
#' @param name sample name
#' @param expermient experiment name
#' 
#' @keywords preprocessing, marker, gene
#' 
#' @return Nothing
do_marker_gene_selection <- function(seu, out_path=NULL){
  markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  if(!is.null(out_path)){ saveRDS(markers, paste0(out_path, ".markers.RDS"))}
  return(markers)
}


##########################################################################


#' do_subsetting
#' Subset samples according to experimental condition
#'
#' @param exp_design Experiment design table in TSV format
#' @param column Column with the condition used for subsetting
#' 
#' @keywords preprocessing, subsetting, integration
#' 
#' @return List of vectors with sample names
do_subsetting <- function(exp_design, column){
  exp_design <- read.csv(exp_design, sep = "\t")
  subsets <- split(exp_design$code, exp_design[[column]])
  return(subsets)
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
  exp_design <- read.csv(exp_design, sep = "\t")
  exp_design <- as.list(exp_design[exp_design$code == name,])
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
merge_seurat <- function(project_name, samples, exp_design,
                                 count_path, suffix=''){
  full_paths <- Sys.glob(paste(count_path, suffix, sep = "/"))
  seu.list <- sapply(samples, function(sample){
    sample_path <- grep(sample, full_paths, value = TRUE)
    d10x <- Read10X(sample_path)
    # Add sample name at the end of cell names
    colnames(d10x) <- paste(sapply(strsplit(colnames(d10x), split="-"),
                                   '[[', 1L), sample, sep="-") 
    seu <- CreateSeuratObject(counts = d10x, project = sample, min.cells = 1,
                              min.features = 1, assay = "scRNAseq")
    seu <- add_exp_design(seu = seu, name = sample, exp_design = exp_design)
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    seu <- subset(seu, subset = nFeature_scRNAseq > 500 & nCount_scRNAseq > 1000
                                & percent.mt < 20)
    return(seu)
    })
  merged_seu <- scCustomize::Merge_Seurat_List(list_seurat = seu.list,
                                               project = project_name)
  return(merged_seu)
}

#' annotate_clusters
#' `annotate_clusters` renames seurat clusters according to dictionary
#'
#' @param seu Non-annotated seu with markers
#' @param anno_table Path to file containing cluster annotation
#' @returns Annotated seu object

annotate_clusters <- function(seu, anno_table) {
  new_clusters <- read.table(file = anno_table, sep = '\t', header = FALSE)[, 1]
  names(new_clusters) <- levels(seu)
  seu[['seurat_annotations']] <- new_clusters
  Idents(seu) <- "seurat_annotations"
  return(seu)
}

#' get_sc_DEGs
#' `get_sc_DEGs` performs differential expression analysis on a seurat object.
#'
#' @param seu Seurat object on which DEG analysis will be performed
#' @param cond A string. Condition by which to perform DEG analysis.
#' @returns A data frame. Rows contain markers, and columns are p_value,
#' average log2FC, pct.1 ? and pct.2 ?. Also adjusted p-value.

get_sc_DEGs <- function(seu, cond) {
  if(!cond %in% names(seu@meta.data)) {
    stop("Specified condition does not exist in seurat object metadata.
          DEG analysis impossible.")
  }
  Idents(seu) <- cond
  conds <- unique(seu[[cond]][[1]])
  if(length(cond) != 2) {
    stop(paste0("Error: Two groups must be supplied for DEG analysis. Provided ",
                length(cond)))
  }
  markers <- FindMarkers(seu, ident.1 = conds[1], ident.2 = conds[2])
  return(markers)
}

#' subset_seurat
#' `subset_seurat` allows subsetting a seurat object without requiring
#' literal strings.
#'
#' @param seu Seurat object on which DEG analysis will be performed
#' @param col Metadata column by which to subset
#' @param value Condition to subset
#' @returns A subsetted seurat object

subset_seurat <- function(seu, col) {
  subset_names <- unique(seu@meta.data[[col]])
  expr <- FetchData(seu, vars = col)
  subs_seu <- list()
  subs_seu <- lapply(subset_names, function (name) {
                                  subs_seu[[name]] <- seu[, which(expr == name)]
    })
  names(subs_seu) <- subset_names
  return(subs_seu)
}

#' do_harmony
#' Perform integration of a merged Seurat object with Harmony
#'
#' @param seu Merged Seurat object
#' @param exp_cond Seurat metadata column with the sample names
#' 
#' @keywords preprocessing, integration
#' 
#' @return Multi-sample integrated Seurat object with Harmony embeddings available
do_harmony <- function(seu, exp_cond){
  seu <- harmony::RunHarmony(seu, exp_cond)
  return(seu)
}


##########################################################################


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
#' @param clusters_annotation Path to clusters annotation file.
#'
#' @keywords preprocessing, main
#' 
#' @return Seurat object
analyze_seurat <- function(raw_seu, out_path = NULL, minqcfeats, percentmt,
                           normalmethod, scalefactor, hvgs, ndims, resolution,
                           dimreds_to_do, embeddings_to_use, integrate = FALSE,
                           clusters_annotation = NULL){
  raw_seu <- do_qc(seu = raw_seu, minqcfeats = minqcfeats, 
                   percentmt = percentmt)
  # Save before version
  if (!is.null(out_path)) saveRDS(raw_seu, paste0(out_path, ".before.seu.RDS")) 
  seu <- subset(raw_seu, subset = QC != 'High_MT,Low_nFeature')
  seu <- NormalizeData(seu, normalization.method = normalmethod,
                       scale.factor = scalefactor, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = hvgs)
  seu <- ScaleData(seu, features = rownames(seu))
  # dimreds_to_do: PCA/UMAP/tSNE
  seu <- do_dimred(seu = seu, ndims = ndims, dimreds = dimreds_to_do)   
  if(integrate) { # Harmony integration and remaining dimreds 
    seu <- do_harmony(seu = seu, exp_cond = "code")
    seu <- do_dimred(seu = seu, ndims = ndims, dimreds = c("tsne", "umap"),
                     reduction = embeddings_to_use)
  }
  seu <- do_clustering(seu = seu, ndims = ndims, resolution = resolution,
                       reduction = embeddings_to_use)
  seu <- JoinLayers(seu)
  markers <- do_marker_gene_selection(seu = seu, out_path = out_path)
  if(!is.null(out_path)) saveRDS(seu, paste0(out_path, ".seu.RDS"))
  return(list(seu = seu, raw_seu = raw_seu, markers = markers))
}

#' integrate_seurat
#' `integrate_seurat` subsets a seurat object and separately analyzes each
#' condition.
#'
#' @param seu seurat object to analyze
#' @param columns experimental design columns by which to integrate
#' returns a list containing the result of each comparison

integrate_seurat <- function(seu, column) {
  values <- unique(column)
  subset <- subset_seurat(seu, column)
  all_subset <- analyze_seurat(subset)
  return(all_subset)
}


##########################################################################


#' write_preprocessing_report
#' Write preprocessing HTML report
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

write_integration_report <- function(seu, template, intermediate_files,
                                     markers_general, markers_specific) {
  return('placeholder')
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

