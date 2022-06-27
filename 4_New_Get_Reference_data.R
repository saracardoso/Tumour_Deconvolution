atlas_file = './Data/CRC_atlas/CRCatlas.h5Seurat'
source('./code/R/EstimateProportions_Algorithms.R')





# ---
# - Load Atlas
# ---

crc_atlas = SeuratDisk::LoadH5Seurat(atlas_file)





# ---
# - Get initial reference dataset
# ---

# 1. Remove cells from normal samples:
crc_atlas = subset(crc_atlas, cells=rownames(crc_atlas[[]])[crc_atlas$Sample.State=='Tumor' & 
                                                              crc_atlas$Annotation_1!='Normal Epithelial cells'])
invisible(gc())

# 2. Add annotation according to ground truth
crc_atlas[['Deconv_cellTypes']] = rep('', dim(crc_atlas[[]])[1])
crc_atlas$Deconv_cellTypes[crc_atlas$Annotation_1%in%c('Mast cells', 'Unknown', 'IgG+ Plasma cells', 'IgA+ Plasma cells',
                                                       'Unconventional Tcells', 'pDCs', 'cDCs')] = 'Other cells'
crc_atlas$Deconv_cellTypes[crc_atlas$Annotation_3%in%c('LTi cells', 'Double-Negative Tcells')] = 'Other cells'
crc_atlas$Deconv_cellTypes[crc_atlas$Annotation_1=='Tumour Epithelial cells'] = 'Cancer cells'
crc_atlas$Deconv_cellTypes[crc_atlas$Annotation_0=='Stromal cells'] = 'Stromal cells'
crc_atlas$Deconv_cellTypes[crc_atlas$Annotation_1%in%c('Naive Bcells', 'Memory Bcells', 'Proliferative Bcells')] = 'Bcells'
crc_atlas$Deconv_cellTypes[crc_atlas$Annotation_1%in%c('Anti-inflammatory macro/mono',
                                                       'Pro-inflammatory macro/mono')] = 'Macro/mono Lineage'
crc_atlas$Deconv_cellTypes[crc_atlas$Annotation_3=='NK cells'] = 'NK cells'
crc_atlas$Deconv_cellTypes[crc_atlas$Annotation_4%in%c('Proliferative CD4 Tcells', 'Proliferative CD8 Tcells')] = 'Proliferative Tcells'
crc_atlas$Deconv_cellTypes[crc_atlas$Annotation_4%in%c('Naive CD4 Tcells', 'IL22+ CD4 Tcells', 'Memory CD4 Tcells', 'IL17+ CD4 Tcells',
                                                       'Follicular CD4 Tcells')] = 'CD4 Tcells'
crc_atlas$Deconv_cellTypes[crc_atlas$Annotation_4=='Regulatory CD4 Tcells'] = 'Regulatory CD4 Tcells'
crc_atlas$Deconv_cellTypes[crc_atlas$Annotation_4%in%c('Naive CD8 Tcells', 'CXCL13+ CD8 Tcells',
                                                       'Memory CD8 Tcells', 'Cytotoxic CD8 Tcells')] = 'CD8 Tcells'

# 3. Subset data. Pool 1 500 cells from each group, except 'other cells'. If that group has less than
#   1 500 cells, then all cells are used.
ct_nCells = table(crc_atlas$Deconv_cellTypes)
cells_use = c()
for(ct in unique(crc_atlas$Deconv_cellTypes)){
  if(!ct%in%'Other cells'){
    n_cells = ct_nCells[ct]
    if(ct_nCells[ct] > 1200) n_cells = 1200
    ct_cells = rownames(crc_atlas[[]])[crc_atlas$Deconv_cellTypes==ct]
    cells_use = c(cells_use, sample(ct_cells, n_cells))
  }
}
deconv_reference = subset(crc_atlas, cells=cells_use)

# 4. Visually check cell-types:
Seurat::DimPlot(deconv_reference, group.by='Deconv_cellTypes')

# 5. Save seurat object
SeuratDisk::SaveH5Seurat(deconv_reference, './Data/CRC_reference/new/CRC_initial.h5Seurat')
remove(crc_atlas)
invisible(gc())





# ---
# - Get Metadata file
# ---

metadata_deconv_reference = deconv_reference@meta.data[, c('Dataset', 'Individual', 'Sample', 'Deconv_cellTypes')]
rownames(metadata_deconv_reference) = gsub('[-]', '.', rownames(metadata_deconv_reference))
write.csv(metadata_deconv_reference, './Data/CRC_reference/new/metadata.csv')




# ---
# - Calculate markers between cell-types
# ---

# 1. Get markers:
Seurat::Idents(deconv_reference) = 'Deconv_cellTypes'
markers = Seurat::FindAllMarkers(deconv_reference, assay='RNA', slot='data',
                                 logfc.threshold=0.8, min.pct=0.3, only.pos=TRUE)

# 2. Convert results into a list. Only store genes whose adjusted p-value is smaller than 0.01:
markers_list = c()
for(ct in unique(markers$cluster)){
  markers_list[[ct]] = markers[markers$cluster==ct & markers$p_val_adj < 0.01, 'gene']
}

# 3. Save list of markers into a json file:
jsonlite::write_json(markers_list, './Data/CRC_reference/new/markers.json')





# ---
# - Get raw counts
# ---

# 1. Raw counts
# 1.1. Get matrix:
raw_counts = Seurat::GetAssayData(deconv_reference, assay='RNA', slot='counts')
colnames(raw_counts) = gsub('-', '.', colnames(raw_counts))
# 1.2. Eliminate genes with 0 counts in all cells
raw_counts = raw_counts[Matrix::rowSums(raw_counts)!=0, ]





# ---
# - Raw counts | Without bias correction 
# ---

nonprolT = rownames(metadata_deconv_reference)[metadata_deconv_reference$Deconv_cellTypes!='Proliferative Tcells']

# 1. Store matrix in csv format:
write.csv(as.matrix(raw_counts), './Data/CRC_reference/new/UB_matrix_CSV.csv')

# 2. Store matrix in CIBERSORTx format (txt):
raw_counts_cbsrtx = as.matrix(raw_counts)
colnames(raw_counts_cbsrtx) =  make.unique(metadata_deconv_reference[,'Deconv_cellTypes'])
colnames(raw_counts_cbsrtx)[grep('[.]', colnames(raw_counts_cbsrtx), invert=T)] =
  paste(colnames(raw_counts_cbsrtx)[grep('[.]', colnames(raw_counts_cbsrtx), invert=T)], '.0', sep='')
colnames(raw_counts_cbsrtx)[grep(' ', colnames(raw_counts_cbsrtx))] =
  gsub(' ', '_', colnames(raw_counts_cbsrtx)[grep(' ', colnames(raw_counts_cbsrtx))])
write.table(raw_counts_cbsrtx, './Data/CRC_reference/new/UB_matrix_TXT.txt', col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')
# Matrix without the proliferative cells:
raw_counts_cbsrtx = as.matrix(raw_counts[,nonprolT])
colnames(raw_counts_cbsrtx) =  make.unique(metadata_deconv_reference[nonprolT,'Deconv_cellTypes'])
colnames(raw_counts_cbsrtx)[grep('[.]', colnames(raw_counts_cbsrtx), invert=T)] =
  paste(colnames(raw_counts_cbsrtx)[grep('[.]', colnames(raw_counts_cbsrtx), invert=T)], '.0', sep='')
colnames(raw_counts_cbsrtx)[grep(' ', colnames(raw_counts_cbsrtx))] =
  gsub(' ', '_', colnames(raw_counts_cbsrtx)[grep(' ', colnames(raw_counts_cbsrtx))])
write.table(raw_counts_cbsrtx, './Data/CRC_reference/new/UB_matrix_woProl_TXT.txt', col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')


# 3. BSeqSC reference matrix:
bseqsc_S = create_reference_matrix_BseqSC(list(data=as.matrix(raw_counts), metadata=metadata_deconv_reference), 'Deconv_cellTypes',
                                          markers=markers_list, ct_scale=T, rna_bias=FALSE, bias_vec=NULL)
colnames(bseqsc_S)[grep(' ', colnames(bseqsc_S))] = gsub(' ', '_', colnames(bseqsc_S)[grep(' ', colnames(bseqsc_S))])
write.table(bseqsc_S, './Data/CRC_reference/new/UB_matrix_BseqSC.txt', col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')
invisible(gc())
# Matrix without the proliferative cells:
bseqsc_S = create_reference_matrix_BseqSC(list(data=as.matrix(raw_counts[, nonprolT]), metadata=metadata_deconv_reference[nonprolT,]), 'Deconv_cellTypes',
                                          markers=markers_list[-7], ct_scale=T, rna_bias=FALSE, bias_vec=NULL)
invisible(gc())
colnames(bseqsc_S)[grep(' ', colnames(bseqsc_S))] = gsub(' ', '_', colnames(bseqsc_S)[grep(' ', colnames(bseqsc_S))])
write.table(bseqsc_S, './Data/CRC_reference/new/UB_matrix_woProl_BseqSC.txt', col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')

# 4. Scaden train file:
# 4.1. With proliferative cells
sc_processed_scaden = processing_scData_scaden(list(data=as.matrix(raw_counts), metadata=metadata_deconv_reference), 'Deconv_cellTypes',
                                               dataset_variable=NULL)
invisible(gc())
create_train_dataFile_scaden(sc_processed_scaden, './Data/CRC_reference/new/UB_Scaden.h5ad', n_samples=1000L, n_cells=100L)
invisible(gc())
# 4.2. Without proliferative cells
raw_counts_woProl = raw_counts[, nonprolT]
meta_woProl = metadata_deconv_reference[nonprolT,]
invisible(gc())
sc_processed_scaden = processing_scData_scaden(list(data=as.matrix(raw_counts_woProl),
                                                    metadata=meta_woProl),
                                               'Deconv_cellTypes', dataset_variable=NULL)
invisible(gc())
create_train_dataFile_scaden(sc_processed_scaden, './Data/CRC_reference/new/UB_Scaden_woProl.h5ad', n_samples=1000L, n_cells=100L)
invisible(gc())





# ---
# - Raw counts | With bias correction 
# ---

nonprolT = rownames(metadata_deconv_reference)[metadata_deconv_reference$Deconv_cellTypes!='Proliferative Tcells']

# 1. Calculate bias:
bias_vector = c()
total_gene_counts_cells = Matrix::colSums(raw_counts)
for(ct in unique(metadata_deconv_reference$Deconv_cellTypes)){
  bias_vector = c(bias_vector, mean(total_gene_counts_cells[metadata_deconv_reference$Deconv_cellTypes==ct]))
}
names(bias_vector) = unique(metadata_deconv_reference$Deconv_cellTypes)
# 1.1. Save bias vector:
write.csv(bias_vector, './Data/CRC_reference/new/bias_values.csv')

# 2. Create biased counts:
biased_counts = as.matrix(raw_counts)
for(ct in unique(metadata_deconv_reference$Deconv_cellTypes)){
  message(ct)
  ct_cells = rownames(metadata_deconv_reference)[metadata_deconv_reference$Deconv_cellTypes==ct]
  biased_counts[,ct_cells] = biased_counts[,ct_cells] * bias_vector[ct]
}

# 3. Store matrix in csv format:
write.csv(as.matrix(biased_counts), './Data/CRC_reference/new/B_matrix_CSV.csv')

# 4. Store matrix in CIBERSORTx format (txt):
biased_counts_cbsrtx = as.matrix(biased_counts)
colnames(biased_counts_cbsrtx) =  make.unique(metadata_deconv_reference[,'Deconv_cellTypes'])
colnames(biased_counts_cbsrtx)[grep('[.]', colnames(biased_counts_cbsrtx), invert=T)] =
  paste(colnames(biased_counts_cbsrtx)[grep('[.]', colnames(biased_counts_cbsrtx), invert=T)], '.0', sep='')
colnames(biased_counts_cbsrtx)[grep(' ', colnames(biased_counts_cbsrtx))] =
  gsub(' ', '_', colnames(biased_counts_cbsrtx)[grep(' ', colnames(biased_counts_cbsrtx))])
write.table(biased_counts_cbsrtx, './Data/CRC_reference/new/B_matrix_TXT.txt', col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')
invisible(gc())
# Make matrix without proliferative Tcells
biased_counts_cbsrtx = as.matrix(biased_counts[, nonprolT])
colnames(biased_counts_cbsrtx) =  make.unique(metadata_deconv_reference[nonprolT, 'Deconv_cellTypes'])
colnames(biased_counts_cbsrtx)[grep('[.]', colnames(biased_counts_cbsrtx), invert=T)] =
  paste(colnames(biased_counts_cbsrtx)[grep('[.]', colnames(biased_counts_cbsrtx), invert=T)], '.0', sep='')
colnames(biased_counts_cbsrtx)[grep(' ', colnames(biased_counts_cbsrtx))] =
  gsub(' ', '_', colnames(biased_counts_cbsrtx)[grep(' ', colnames(biased_counts_cbsrtx))])
write.table(biased_counts_cbsrtx, './Data/CRC_reference/new/B_matrix_woProl_TXT.txt', col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')
invisible(gc())

# 5. BSeqSC reference matrix:
bseqsc_SB = create_reference_matrix_BseqSC(list(data=as.matrix(raw_counts), metadata=metadata_deconv_reference), 'Deconv_cellTypes',
                                           markers=markers_list, ct_scale=T, rna_bias=TRUE, bias_vec=bias_vector)
invisible(gc())
colnames(bseqsc_SB)[grep(' ', colnames(bseqsc_SB))] = gsub(' ', '_', colnames(bseqsc_SB)[grep(' ', colnames(bseqsc_SB))])
write.table(bseqsc_SB, './Data/CRC_reference/new/B_matrix_BseqSC.txt', col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')
invisible(gc())
# Make matrix without prolifrative Tcells
bseqsc_SB = create_reference_matrix_BseqSC(list(data=as.matrix(raw_counts[, nonprolT]), metadata=metadata_deconv_reference[nonprolT,]), 'Deconv_cellTypes',
                                           markers=markers_list[-7], ct_scale=T, rna_bias=TRUE, bias_vec=bias_vector)
invisible(gc())
colnames(bseqsc_SB)[grep(' ', colnames(bseqsc_SB))] = gsub(' ', '_', colnames(bseqsc_SB)[grep(' ', colnames(bseqsc_SB))])
write.table(bseqsc_SB, './Data/CRC_reference/new/B_matrix_woProl_BseqSC.txt', col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')
invisible(gc())

# 6. Scaden train file:
# 6.1. With proliferative cells
sc_processed_scaden = processing_scData_scaden(list(data=as.matrix(biased_counts), metadata=metadata_deconv_reference), 'Deconv_cellTypes',
                                               dataset_variable=NULL)
invisible(gc())
create_train_dataFile_scaden(sc_processed_scaden, './Data/CRC_reference/new/B_Scaden.h5ad', n_samples=1000L, n_cells=100L)
invisible(gc())
# 6.2. Without proliferative cells
biased_counts_woProl = as.matrix(Seurat::GetAssayData(subset(deconv_reference, cells=nonprolT),
                                                      assay='RNA', slot='counts'))
invisible(gc())
for(ct in unique(metadata_deconv_reference$Deconv_cellTypes)){
  message(ct)
  if(ct=='Proliferatice Tcells') next
  ct_cells = rownames(metadata_deconv_reference)[metadata_deconv_reference$Deconv_cellTypes==ct]
  biased_counts_woProl[,ct_cells] = biased_counts_woProl[,ct_cells] * bias_vector[ct]
}
invisible(gc())
sc_processed_scaden = processing_scData_scaden(list(data=as.matrix(raw_counts_woProl),
                                                    metadata=metadata_deconv_reference[nonprolT,]),
                                               'Deconv_cellTypes', dataset_variable=NULL)
invisible(gc())
create_train_dataFile_scaden(sc_processed_scaden, './Data/CRC_reference/new/B_Scaden_woProl.h5ad', n_samples=1000L, n_cells=100L)
invisible(gc())


















