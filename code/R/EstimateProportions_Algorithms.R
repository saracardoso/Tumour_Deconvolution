dwls_modified_file = './code/R/DWLS_modified.R'
ddlsorter_modified_file = './code/R/DigitalDLSorter_simSCProfiles_modified.R'
rna_correction_file = './code/R/rna_correction.R'
scaden_bulk_simulation_modified_file = './code/inst/SCADEN_bulk_simulation_modified.py'
scaden_create_h5adFile_modified_file = './code/inst/SCADEN_create_h5ad_file_modified.py'
source(dwls_modified_file)
source(ddlsorter_modified_file)
source(rna_correction_file)





# -------------------------------------------
# --- CREATE REFERENCE MATRIX FOR Bseq-SC ---
# -------------------------------------------

create_reference_matrix_BseqSC = function(sc.reference, cellType.var, markers=NULL, ct_scale=T, rna_bias=FALSE, bias_vec=NULL){
  scReference = prepare_expressionSet_scReferences(sc.reference, cellType.var)
  scReference_matrix = bseqsc::bseqsc_basis(scReference, markers, 'cellType', 'sampleID', ct_scale)
  if(rna_bias){
    if(is.null(bias_vec)) stop('bias_vec argument missing')
      scReference_matrix = correct_reference_mRNABias(scReference_matrix, correction_vector=bias_vec)
    }
  return(scReference_matrix)
}





# --------------------------------------
# --- PREPARE DATA FOR DECONVOLUTION ---
# --------------------------------------


#---
#- Expression Set
#---

prepare_expressionSet_scReferences = function(sc.reference, cellType.var){
  # - phenoData metadata of scRNAseq:
  pheno_metadata = data.frame(labelDescription= c("Sample ID", "Subject Name", "Cell Type ID", "Cell Type Name"),
                              row.names=c("sampleID", "SubjectName", "cellTypeID", "cellType"))
  # - phenoData data and data matrix of scRNAseq:
  pheno_data =  data.frame(sampleID=sc.reference$metadata$Sample, SubjectName=sc.reference$metadata$Sample,
                           cellTypeID=rownames(sc.reference$metadata), cellType=sc.reference$metadata[,cellType.var],
                           row.names=rownames(sc.reference$metadata))
  scData_matrix = sc.reference$data
  colnames(scData_matrix) = row.names(pheno_data)
  scReference = Biobase::ExpressionSet(data.matrix(scData_matrix), phenoData = Biobase::AnnotatedDataFrame(pheno_data, pheno_metadata))
  
  return(scReference)
}


#---
#- DigitalDLSorter object
#---

prepare_singlecellexperiment_scReferences = function(sc.reference, cellType.var, min_cells = 0, min_counts = 0, project_name = 'DDLS'){
  sce = SingleCellExperiment::SingleCellExperiment(as.matrix(sc.reference$data),
                                                   colData = data.frame(Cell_ID = rownames(sc.reference$metadata),
                                                                        Cell_Type = sc.reference$metadata[,cellType.var]),
                                                   rowData = data.frame(Gene_ID = rownames(sc.reference$data)))
  
  DDLSorter_obj = digitalDLSorteR::loadSCProfiles(single.cell = sce, cell.ID.column = "Cell_ID", gene.ID.column = "Gene_ID",
                                                  min.cells = min_cells, min.counts = min_counts, project = project_name)

  return(DDLSorter_obj)
}






# ---------------------
# --- DECONVOLUTION ---
# ---------------------


#---
#- DWLS
#---

deconv_dwls = function(bulk_matrix, scReference, cell_type_variable, diff_cutoff=0.5, pval_cutoff=0.01, #de_method='seurat',
                       rna_bias=FALSE, bias_vec=NULL){
  
  message('Building Signature Matrix\n')
  #if(de_method=='seurat')
  dwls_signature_matrix = buildSignatureMatrixUsingSeurat(scdata=scReference$data, id=scReference$metadata[,cell_type_variable],
                                                          diff.cutoff=diff_cutoff, pval.cutoff=pval_cutoff)
  #else if(de_method=='mast')
  #  dwls_signature_matrix = buildSignatureMatrixMAST(scdata=scReference$data, id=scReference$metadata[,cell_type_variable],
  #                                                   diff.cutoff=diff_cutoff, pval.cutoff=pval_cutoff)
  #else stop('Invalid differential expression method (de_expression). Those available are: seurat.')# and mast.')
  DWLS_proportions = c()
  cat('Calculating Proportions\n')
  for(sample_name in colnames(bulk_matrix)){
    cat('- sample', sample_name, '\n')
    sample_data = bulk_matrix[,sample_name]
    names(sample_data) = rownames(bulk_matrix)
    if(rna_bias){
      if(is.null(bias_vec)) stop('bias_vec argument missing')
      dwls_signature_matrix = correct_reference_mRNABias(dwls_signature_matrix, correction_vector=bias_vec)
    } 
    dwls_data = trimData(dwls_signature_matrix, sample_data)
    DWLS_proportions = cbind(DWLS_proportions, solveDampenedWLS(dwls_data$sig, dwls_data$bulk))
  }
  rownames(DWLS_proportions) = colnames(dwls_signature_matrix)
  colnames(DWLS_proportions) = colnames(bulk_matrix)
  
  cat('Done!\n')
  return(DWLS_proportions)
}



#---
#- MOMF
#---

deconv_momf = function(bulk_matrix, scReference, cell_type_variable, method='KL', rho=2, num_iter=5000,
                       rna_bias=FALSE, bias_vec=NULL){
  
  message('Filtering for common genes between bulk and reference data')
  keep_genes = intersect(rownames(scReference$data), rownames(bulk_matrix))
  scReference$data = scReference$data[keep_genes,]
  bulk_matrix = bulk_matrix[keep_genes,]
  
  message('Building Signature Matrix')
  momf_signature_matrix = MOMF::momf.computeRef(scReference$data, scReference$metadata[,cell_type_variable])
  momf_GL = list(X1=t(scReference$data), X2=t(bulk_matrix))
  
  message('Calculating Proportions')
  if(rna_bias){
    if(is.null(bias_vec)) stop('bias_vec argument missing')
    momf_signature_matrix = correct_reference_mRNABias(momf_signature_matrix, correction_vector=bias_vec)
  }
  MOMF_proportions = MOMF::momf.fit(momf_GL, DataPriorU=momf_signature_matrix, method=method, rho=rho, num_iter=num_iter)
  
  message('Done!')
  return(t(MOMF_proportions$cell.prop))
}



#---
#- CPM
#---

deconv_cpm = function(bulk_matrix, scReference, cell_type_variable, cell_space=NULL,
                      neighborhood_size=10, min_selection=5){
  
  cat('Calculating Proportions\n')
  avlbl_cores = future::availableCores()
  if(avlbl_cores>1) n_cores = avlbl_cores -1
  else n_cores = 1
  CPM_proportions = scBio::CPM(scReference$data, scReference$metadata[,cell_type_variable], bulk_matrix, as.matrix(cell_space),
                               no_cores=n_cores, quantifyTypes=TRUE, typeTransformation=TRUE,
                               neighborhoodSize=neighborhood_size, minSelection=min_selection)
  
  cat('Done!\n')
  return(t(CPM_proportions$cellTypePredictions))
}



#---
#- AutoGeneS
#---

# model can be: nusvr, nnls, linear
# weights can be: minimum correlation (-1,0), minimum correlation and maximum distance: (-1,1), maximum distance: (0,1)
deconv_autogenes = function(bulk_matrix, scReference, cell_type_variable, model='nnls', gene_markers=NULL, mode='fixed',
                            n_iterations=5000L, n_genes=400L, weights = list(-1,0)
                            ){
  # Get necessary python modules:
  reticulate::py_config()
  autogenes_module = reticulate::import('autogenes')
  anndata_module = reticulate::import('anndata')

  # Initialize deconvolution object:
  scData = anndata_module$AnnData(scReference$data)
  scData = scData$transpose()
  scData$obs[, 'celltype'] = scReference$metadata[[cell_type_variable]]
  if(!is.null(gene_markers)){
    higly_variable_boolean = rep(FALSE, dim(scReference$data)[1])
    genes_to_keep = na.omit(match(unique(unlist(gene_markers)), rownames(scReference$data)))
    higly_variable_boolean[genes_to_keep] = TRUE
    scData$var[, 'highly_variable'] = higly_variable_boolean
    autogenes_module$init(scData, use_highly_variable=TRUE, celltype_key='celltype')
  }
  else autogenes_module$init(scData, use_highly_variable=FALSE, celltype_key='celltype')
  
  bData = anndata_module$AnnData(bulk_matrix)
  bData = bData$transpose()
  
  # Get the best set of n_genes genes with lowest correlation between cell types:
  cat('Building Signature Matrix\n')
  if(mode=='fixed') autogenes_module$optimize(ngen=n_iterations, nfeatures=n_genes, seed=0L, mode=mode)
  else autogenes_module$optimize(ngen=n_iterations, seed=0L, mode=mode)
  autogenes_module$select(weights=weights)
  
  # Deconvolve bulk data with selected genes:
  cat('Calculating Proportions\n')
  AutoGeneS_proportions = autogenes_module$deconvolve(bData, model=model)
  
  # Normalize proportions (i.e., make proportions in each sample sum to 1)
  AutoGeneS_proportions_norm = AutoGeneS_proportions
  for(row in 1:dim(AutoGeneS_proportions)[1]){
    AutoGeneS_proportions_norm[row, ] = AutoGeneS_proportions[row, ] / sum(AutoGeneS_proportions[row, ])
  }
  rownames(AutoGeneS_proportions_norm) = colnames(bulk_matrix) # Samples
  colnames(AutoGeneS_proportions_norm) = rownames(autogenes_module$adata()$obs) # Cell types
  
  cat('Done!\n')
  return(as.data.frame(t(AutoGeneS_proportions_norm)))
}



#---
#- MuSiC
#---

# --- MuSiC without pre-grouping of cell-types:

# method ca be: Est.prop.weighted, Est.prop.allgene
deconv_music_woGrouping = function(bulk_data_expressionSet, scRNAseq_expressionSet, method='Est.prop.weighted',
                                   gene_markers=NULL, n_iterations=1000, center=FALSE, normalize=FALSE, bias_vec=NULL){
  library(Biobase)
  library(xbioc)
  if(is.null(bias_vec)) stop('bias_vec argument missing')
  cell_size = data.frame(x=names(bias_vec), y=bias_vec)
  music_proportions = MuSiC::music_prop(bulk.eset=bulk_data_expressionSet, sc.eset=scRNAseq_expressionSet, markers=gene_markers,
                                        clusters='cellType', samples='sampleID', verbose=T,
                                        iter.max=n_iterations, centered=center, normalize=normalize, cell_size=cell_size)
  return(t(music_proportions[[method]]))
}


# --- MuSic with pre-grouping of cell-types:

pre_group_cell_types_music = function(scRef, cell_type_variable, remove_zero=TRUE, gene_markers=NULL,
                                      distance_method='euclidean', agglomeration_method='complete', cut_off=c('design_matrix', 'abundance'),
                                      rna_bias=FALSE, bias_vec=NULL){
  library(Biobase)
  library(xbioc)
  # calculation of groups for MuSiC's pre-grouping algorithm
  scReference = prepare_expressionSet_scReferences(scRef, cell_type_variable)
  if(rna_bias){
    if(is.null(bias_vec)) stop('bias_vec argument missing')
    cell_size = data.frame(x=names(bias_vec), y=bias_vec)
  }
  else cell_size = NULL
  design_matrix = MuSiC::music_basis(scReference, non.zero=remove_zero, markers=gene_markers,
                                     clusters='cellType', samples='sampleID', cell_size=cell_size)
  # Choose matrix to cut off:
  if(length(cut_off)>2) stop('Invalid cut_off argument.')
  else if(length(cut_off)==2) cut_off = 'abundance'
  else if(cut_off%in%c('design_matrix', 'abundance')) cut_off = cut_off
  else stop('Invalid cut_off argument.')
  
  # Calculate clusters:
  if(cut_off == 'design_matrix'){
    # Hierarchical clustering of design matrix:
    d = dist(t(log(design_matrix$Disgn.mtx + 1e-6)), method=distance_method)
    hc = hclust(d, method=agglomeration_method)
    title = 'Cluster log(Design Matrix)'
  }
  else{
    # Hierarchical clustering of average relative abundance:
    d = dist(t(log(design_matrix$M.theta + 1e-8)), method=distance_method)
    hc = hclust(d, method=agglomeration_method)
    title = 'Cluster log(Average of Relative Abundance)'
  }

  # Plot:
  plot(hc, cex=0.6, hang=-1, main=title)
  # Create clusters:
  threshold = summary(d)[5]
  y=rect.hclust(hc, h=threshold)
  clusters = list()
  for(i in 1:length(y)){
    clusters[[paste('C', i, sep='')]] = names(y[[i]])
  }
  return(clusters)
}

calculate_clusters_markers = function(scRNAseq_expressionSet, celltypes_cluster_list,
                                      min_log2FC=0.25, pval=0.01, only_pos=FALSE){
  #Create seurat object:
  exprObj = Seurat::CreateSeuratObject(as.data.frame(exprs(scRNAseq_expressionSet)), project = "DE")
  exprObj = Seurat::SetIdent(exprObj,value=as.vector(scRNAseq_expressionSet$cellType))
  # Normalize data:
  exprObj = Seurat::SCTransform(exprObj)
  # For each cluster where more than one cell-type is present, calculate gene markers between these cell-types:
  res=list()
  for(cl in names(celltypes_cluster_list)){
    message(paste('- Cluster:', cl, sep=' '))
    if(length(celltypes_cluster_list[[cl]])>1){
      markers = c()
      for(cell_type in celltypes_cluster_list[[cl]]){
        message('--', cell_type, ' ')
        other_cells = celltypes_cluster_list[[cl]][-match(cell_type, celltypes_cluster_list[[cl]])]
        cell_type_markers = Seurat::FindMarkers(exprObj, ident.1=cell_type, ident.2=other_cells, slot='counts',
                                                test.use='wilcox', logfc.threshold=min_log2FC, min.pct=0.3, only.pos=only_pos)
        markers = c(markers, unique(rownames(cell_type_markers)[cell_type_markers$p_val_adj<pval]))
      }
      res[[cl]] = markers
    }
    else{
      res[cl] = list(NULL)
      message('Only one cell-type in this cluster')
    }
  }
  return(res)
}

deconv_music_wGrouping = function(bulk_data_expressionSet, scRNAseq_expressionSet, celltypes_cluster_list=NULL,
                                  markers_test='wilcox', markers_min_log2FC=0.25, markers_pval=0.01, markers_only_pos=FALSE,
                                  n_iterations=1000, center=FALSE, normalize=FALSE){
  library(Biobase)
  library(xbioc)
  if(is.null(celltypes_cluster_list)) stop('celltypes_cluster_list missing.')
  if(sum(is.na(match(unlist(celltypes_cluster_list), scRNAseq_expressionSet$cellType)))>0)
    stop('There are cell-types missing in the celltypes_cluster_list.')
  
  
  # Set clusterType variable in sc expressionSet:
  clusterType = as.character(scRNAseq_expressionSet$cellType)
  for(cl in 1:length(celltypes_cluster_list)){
    clusterType[clusterType%in%celltypes_cluster_list[[cl]]] = names(celltypes_cluster_list)[cl]
  }
  pData(scRNAseq_expressionSet)$clusterType = factor(clusterType, levels=names(celltypes_cluster_list))
  
  # Calculate markers for each cluster:
  message('Calculating gene markers for each cluster with more than one cell-type:\n')
  group_markers = calculate_clusters_markers(scRNAseq_expressionSet, celltypes_cluster_list,
                                             markers_min_log2FC, markers_pval, markers_only_pos)
  
  # Calculate Proportions:
  message('Calculating Proportions\n')
  music_proportions = MuSiC::music_prop.cluster(bulk.eset=bulk_data_expressionSet, sc.eset=scRNAseq_expressionSet,
                                                group.markers=group_markers, groups='clusterType', clusters='cellType', samples='sampleID',
                                                clusters.type=celltypes_cluster_list, verbose=T, iter.max=n_iterations,
                                                centered=center, normalize=normalize)
  return(t(music_proportions$Est.prop.weighted.cluster))
}



#---
#- SCDC
#---

# --- SCDC with only one sc reference dataset:

deconv_scdc_single_ref = function(bulk_data_expressionSet, scRNAseq_expressionSet, n_iterations=1000){
  SCDC_proportions = SCDC::SCDC_prop(bulk_data_expressionSet, scRNAseq_expressionSet, 'cellType', 'sampleID',
                                     unique(scRNAseq_expressionSet$cellType), iter.max=n_iterations)
  return(t(SCDC_proportions$prop.est.mvw))
}


# --- SCDC with multiple sc reference datasets:

deconv_scdc_multiple_refs = function(bulk_data_expressionSet, scRNAseq_expressionSet_list, n_iterations=1000,
                                     use_gridSearch=FALSE, search_length=0.05, weights_method=NULL){
  if(use_gridSearch) available_weight_methods = c('inverse SSE','inverse SAE','inverse RMSD','LAD','NNLS','mAD_Y','RMSD_Y','Spearman_Y')
  else available_weight_methods = c('inverse SSE','inverse SAE','inverse RMSD','LAD','NNLS')
  if(!is.null(weights_method))
    if(!weights_method%in%available_weight_methods)
      stop(paste('Invalid weights methods. Available methods:', paste(available_weight_methods, collapse=', '), sep=' '))
    else available_weight_methods = weights_method
  
  cat('Running SCDC for multiple sc References\n')
  SCDC_results = SCDC::SCDC_ENSEMBLE(bulk_data_expressionSet, scRNAseq_expressionSet_list, 'cellType', 'sampleID',
                                         unique(scRNAseq_expressionSet_list[[1]]$cellType), grid.search=use_gridSearch,
                                         search.length=search_length, iter.max=n_iterations)
  
  cat('Calculating the proportions based on the weight method(s)\n')
  n_scLists = length(scRNAseq_expressionSet_list)
  SCDC_proportions = list()
  if(length(available_weight_methods)==1){
    SCDC_proportions = t(SCDC::wt_prop(SCDC_results$w_table[available_weight_methods, 1:n_scLists], SCDC_results$prop.only))
  }
  else{
    for(weight_met in available_weight_methods){
      SCDC_proportions[[weight_met]] = t(SCDC::wt_prop(SCDC_results$w_table[weight_met, 1:n_scLists], SCDC_results$prop.only))
    }
  }
  
  
  return(SCDC_proportions)
}



#---
#- BisqueRNA
#---

deconv_bisquerna = function(bulk_data_expressionSet, scRNAseq_expressionSet, gene_markers=NULL, use_overlap=FALSE, old_cpm=TRUE){
  BisqueRNA_proportions = BisqueRNA::ReferenceBasedDecomposition(bulk_data_expressionSet, scRNAseq_expressionSet,
                                                                 markers=gene_markers, 'cellType', 'sampleID',
                                                                 use.overlap=use_overlap, old.cpm=old_cpm)
  return(BisqueRNA_proportions$bulk.props)
}



#---
#- Scaden
#---

processing_scData_scaden = function(scRef, cell_type_variable, dataset_variable=NULL){
  #Create seurat object:
  exprObj = Seurat::CreateSeuratObject(as.data.frame(scRef$data), project = "DE")
  exprObj = Seurat::SetIdent(exprObj,value=as.vector(scRef$metadata[,cell_type_variable]))
  # Normalize data:
  exprObj = Seurat::NormalizeData(exprObj, normalization.method='RC')
  # Normalized count data:
  norm_count_data = t(as.matrix(Seurat::GetAssayData(exprObj, assay='RNA', slot='data')))
  # Cell type labels:
  cell_type_labels = cbind(rownames(exprObj@meta.data), scRef$metadata[rownames(exprObj@meta.data), cell_type_variable])
  if(!is.null(dataset_variable)) cell_type_labels = cbind(cell_type_labels, exprObj@meta.data[, dataset_variable])
  else cell_type_labels = cbind(cell_type_labels, rep('ds1', dim(cell_type_labels)[1]))
  rownames(cell_type_labels) = rownames(exprObj@meta.data)
  colnames(cell_type_labels) = c('sample_names', 'Celltype', 'ds')
  cell_type_labels = cell_type_labels[rownames(norm_count_data), ]
  
  return(list(data=norm_count_data, cell_types=cell_type_labels))
}

create_train_dataFile_scaden = function(processed_data_scaden, trainData_filename,
                                        n_samples=8000L, n_cells=100L, unknown_cell_types=c('Unknown', 'unknown')){
  n_samples2 = floor(n_samples/2)
  
  # Original Data:
  xs = list()
  ys = list()
  datasets = unique(processed_data_scaden$cell_types[,'ds'])
  message(paste('Datasets:', paste(datasets, collapse=', '), sep=' '))
  for(ds in datasets){
    xs[[ds]] = processed_data_scaden$data[processed_data_scaden$cell_types[,'ds']==ds,]
    ys[[ds]] = processed_data_scaden$cell_types[processed_data_scaden$cell_types[,'ds']==ds,]
  }
  
  # Merge unknown cell-types:
  cell_types = unique(processed_data_scaden$cell_types[,'Celltype'])
  which_unknown = na.omit(match(cell_types, unknown_cell_types))
  if(length(which_unknown)>0){
    for(ds in names(ys)){
      unknown_cells = ys[[ds]][, 'Celltype']%in%unknown_cell_types[which_unknown]
      ys[[ds]][, 'Celltype'][unknown_cells] = 'Unknown'
    }
    cell_types = c(cell_types, 'Unknown')
  }
  message(paste('Available celltypes:', paste(cell_types, collapse=', '), sep=' '))
  
  # Create train datasets:
  tmpx = list()
  tmpy = list()
  create_pseuso_bulk_scaden_code = reticulate::py_run_file(scaden_bulk_simulation_modified_file)
  for(ds in names(xs)){
    message(paste('Subsampling ', ds, '...', sep=''))
    tmp = create_pseuso_bulk_scaden_code$create_subsample_dataset(as.data.frame(xs[[ds]]), as.data.frame(ys[[ds]]),
                                                                  as.integer(n_cells), as.list(cell_types), as.integer(n_samples2))
    tmpx[[ds]] = tmp[[1]]
    tmpy[[ds]] = tmp[[2]]
  }
  
  # Create h5ad file:
  message('Creating h5ad file')
  create_h5adFile_scaden_code = reticulate::py_run_file(scaden_create_h5adFile_modified_file)
  create_h5adFile_scaden_code$create_AnnData_dataset(tmpx, tmpy, cell_types, unknown_cell_types, trainData_filename)
  
}

deconv_scaden = function(bulk_filename, reference_h5ad_filename, cell_type_variable, dir_results=tempdir(),
                         datasets='', min_expression=0.1, batch_size=128L, learning_rate=0.0001, n_steps=1000L){
  # Get necessary python modules:
  reticulate::py_config()
  scaden_module = reticulate::import('scaden')
  anndata_module = reticulate::import('anndata')
  
  # Processing:
  cat('Filtering bulk and sc reference(s) for common genes\n')
  processed_reference_filename = paste(dir_results, 'processed_reference.h5ad', sep='/')
  scaden_module$scaden_main$processing(bulk_filename, reference_h5ad_filename, processed_reference_filename, min_expression)
  
  # Training:
  cat('Training sc reference(s)\n')
  scaden_module$scaden_main$training(processed_reference_filename, datasets, dir_results, batch_size, learning_rate, n_steps)
  
  # Prediction:
  cat('Predicting bulk data proportions\n')
  predictions_filename = paste(dir_results, 'estimated_proportions.txt', sep='/')
  scaden_module$scaden_main$prediction(dir_results, bulk_filename, predictions_filename)
  
  # Read predictions file and return it as a result:
  proportions_scaden = read.table(predictions_filename, sep='\t', header=TRUE, row.names=1)
  colnames(proportions_scaden) = gsub('[.]', ' ', colnames(proportions_scaden))
  
  return(t(proportions_scaden))
}



#---
#- DigitalDLSorter
#---

# add_method: fixed or limit. If fixed, all cell-types in cell.types will be added ncells. If limit, all cell-types in cell.types will be
# added cells to reach the limit of ncells per cell-type
simulate_cells = function(DDLSorter_obj, ncells, cell.types, subset.cells=NULL, proportional=TRUE, threads=1, add_method='fixed'){
  # Estimate parameters:
  DDLSorter_objNew = digitalDLSorteR::estimateZinbwaveParams(object = DDLSorter_obj, cell.ID.column = "Cell_ID", gene.ID.column = "Gene_ID",
                                                             cell.type.column = "Cell_Type", subset.cells = subset.cells,
                                                             proportional=proportional, threads = threads)
  # Simulate cells with estimated parameters:
  if(add_method=='fixed'){
    DDLSorter_objNew = digitalDLSorteR::simSCProfiles(object = DDLSorter_objNew, cell.ID.column = "Cell_ID", cell.type.column = "Cell_Type",
                                                      n.cells = ncells, cell.types = cell.types)
  }
  else if(add_method=='limited'){
    ncells_add = ncells - table(DDLSorter_objNew@single.cell.real$Cell_Type)[cell.types]
    DDLSorter_objNew = simSCProfiles_modified(object = DDLSorter_objNew, cell.ID.column = "Cell_ID", cell.type.column = "Cell_Type",
                                              n.cells = ncells_add, cell.types = cell.types)
  }
  return(DDLSorter_objNew)
}

generate_bulk = function(DDLSorter_obj, minmaxprobs,  nsamples = 250, ncells = 100, balanced_cellTypes = FALSE, threads=1){
  # Generate proportions matrix:
  DDLSorter_objNew = digitalDLSorteR::generateBulkCellMatrix(object = DDLSorter_obj, cell.ID.column = "Cell_ID", cell.type.column = "Cell_Type",
                                            prob.design = minmaxprobs, num.bulk.samples = nsamples, n.cells = ncells,
                                            balanced.type.cells=balanced_cellTypes)
  # Generate bulk profile:
  DDLSorter_objNew = digitalDLSorteR::simBulkProfiles(object = DDLSorter_objNew, threads=threads)
  
  return(DDLSorter_objNew)
}

train_DDLSorter_model = function(DDLSorter_obj, combine = "both", batch.size = 64, nepochs = 10, threads = 1, view.metrics.plot = TRUE){
  DDLSorter_objTrained = digitalDLSorteR::trainDigitalDLSorterModel(object = DDLSorter_obj, combine = "both", batch.size = batch.size,
                                                                    num.epochs = nepochs, threads = threads,
                                                                    view.metrics.plot = view.metrics.plot)
  return(DDLSorter_objTrained)
}

deconv_DigitalDLSorter = function(DDLSorter_obj, bulk.data, batch.size = 128, normalize = TRUE, data_name='deconv'){
  # Create summarized experiment for bulk data:
  seBulk = SummarizedExperiment::SummarizedExperiment(assay = list(counts = bulk.data$data))
  # Load bulk data into DDLSorter object:
  DDLSorter_objDeconv = digitalDLSorteR::loadDeconvData(object = DDLSorter_obj, data = seBulk, name.data = data_name)
  # Deconvolute bulk data:
  DDLSorter_objDeconv = digitalDLSorteR::deconvDigitalDLSorterObj(object = DDLSorter_objDeconv, batch.size = batch.size,
                                                                  normalize = normalize, name.data = data_name, verbose = TRUE)
  # Get deconvoluted proportions:
  props_results = t(DDLSorter_objDeconv@deconv.results$deconv)
  return(props_results)
}
