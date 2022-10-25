code_dir = './code/R'
for(file in list.files(code_dir, full.names=TRUE)) source(file)

methods = c('AutoGeneS', 'BisqueRNA', 'BSeqSC', 'CIBERSORTx', 'DigitalDLSorter', 'DWLS', 'MOMF', 'MuSiC_woGrouping', 'MuSiC_wGrouping', 'Scaden', 'SCDC')
cell_types = c('Cancer cells', 'Stromal cells', 'Macro/mono Lineage', 'Bcells', 'CD4 Tcells',
               'Regulatory CD4 Tcells', 'CD8 Tcells', 'Proliferative Tcells', 'NK cells')
cell_types_wOther = c(cell_types, 'Other cells')





# ---
# - Read bulk and reference files
# ---

# 1. Read bulk data:
CRC_bulk = read_bulk_data('./Data/bulk/CRC/NICs/data/data.csv',
                          ground_truth_file='./Data/bulk/CRC/NICs/data/cell_proportions.csv')
names(CRC_bulk) = c('data', 'metadata')
CRC_bulk_file_scaden = './Data/bulk/CRC/NICs/data/data.txt'


# 2. Read reference files:
# 2.1. Reference counts:
references = list()
references$uncorrected = read_scRNAseq_reference('./Data/CRC_reference/UB_matrix_CSV.csv',
                                                 './Data/CRC_reference/metadata.csv')
invisible(gc())
references$corrected = read_scRNAseq_reference('./Data/CRC_reference/B_matrix_CSV.csv',
                                               './Data/CRC_reference/metadata.csv')
invisible(gc())

# 2.2. Gene Markers:
markers =  jsonlite::read_json('./Data/CRC_reference/markers.json', simplifyVector=T)

# 2.3. Scaden train file:
scaden_train_files = list()
scaden_train_files$corrected = list()
scaden_train_files$corrected$wProl = './Data/CRC_reference/B_Scaden.h5ad'
scaden_train_files$corrected$woProl = './Data/CRC_reference/B_Scaden_woProl.h5ad'
scaden_train_files$uncorrected = list()
scaden_train_files$uncorrected$wProl = './Data/CRC_reference/UB_Scaden.h5ad'
scaden_train_files$uncorrected$woProl = './Data/CRC_reference/UB_Scaden_woProl.h5ad'

# 2.4. RNA bias correction vector
rna_correction_df = read.csv('./Data/CRC_reference/bias_values.csv', row.names=1)
rna_correction_vec = rna_correction_df[,1]
names(rna_correction_vec) = rownames(rna_correction_df)





# ---
# - Run methods
# ---

notprol_Tcells = rownames(references$uncorrected$metadata)[references$uncorrected$metadata$Deconv_cellTypes!='Proliferative Tcells']

# 1. AutoGeneS

# 1.1. Run autogenes, with proliferative Tcells
# 1.1.1. autogenes_linear
res = Estimate.Proportions(CRC_bulk, references$uncorrected, 'Deconv_cellTypes', 'AutoGeneS',
                           model='linear', gene.markers=markers)
invisible(gc())
write.csv(res, './Results/proportions/AutoGeneS_linear.csv')
# 1.1.2. autogenes_nnls
res = Estimate.Proportions(CRC_bulk, references$uncorrected, 'Deconv_cellTypes', 'AutoGeneS',
                           model='nnls', gene.markers=markers)
invisible(gc())
write.csv(res, './Results/proportions/AutoGeneS_nnls.csv')
# 1.1.3. autogenes_nusvr
res = Estimate.Proportions(CRC_bulk, references$uncorrected, 'Deconv_cellTypes', 'AutoGeneS',
                           model='nusvr', gene.markers=markers)
invisible(gc())
write.csv(res, './Results/proportions/AutoGeneS_nusvr.csv')

# 1.2. Run autogenes, without proliferative Tcells
# 1.2.1. autogenes_linear
res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                      metadata=references$uncorrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'AutoGeneS', model='linear', gene.markers=markers)
invisible(gc())
write.csv(res, './Results/proportions/woProl_AutoGeneS_linear.csv')
# 1.2.2. autogenes_nnls
res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                      metadata=references$uncorrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'AutoGeneS', model='nnls', gene.markers=markers)
invisible(gc())
write.csv(res, './Results/proportions/woProl_AutoGeneS_nnls.csv')
# 1.2.3. autogenes_nusvr
res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                      metadata=references$uncorrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'AutoGeneS', model='nusvr', gene.markers=markers)
invisible(gc())
write.csv(res, './Results/proportions/woProl_AutoGeneS_nusvr.csv')



# 2. BisqueRNA

# 2.1. Run BisqueRNA, with proliferative Tcells
# 2.1.1. Without correction
res = Estimate.Proportions(CRC_bulk, references$uncorrected, 'Deconv_cellTypes', 'BisqueRNA',
                                gene.markers=markers, use.overlap=FALSE, old.cpm=TRUE)
invisible(gc())
write.csv(res, './Results/proportions/uncorrected/BisqueRNA.csv')
# 2.1.2. With correction before
res = Estimate.Proportions(CRC_bulk, references$corrected, 'Deconv_cellTypes', 'BisqueRNA',
                                gene.markers=markers, use.overlap=FALSE, old.cpm=TRUE)
invisible(gc())
write.csv(res, './Results/proportions/corrected_before/BisqueRNA.csv')
# 2.1.3. With correction after
res = Estimate.Proportions(CRC_bulk, references$uncorrected, 'Deconv_cellTypes', 'BisqueRNA',
                                gene.markers=markers, use.overlap=FALSE, old.cpm=TRUE)
invisible(gc())
res = correct_fractions_mRNABias(res, rna_correction_vec)
write.csv(res, './Results/proportions/corrected_after/BisqueRNA.csv')

# 2.2. Run BisqueRNA, without proliferative Tcells
# 2.2.1. Without correction
res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                          metadata=references$uncorrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'BisqueRNA', gene.markers=markers, use.overlap=FALSE,
                           old.cpm=TRUE)
invisible(gc())
write.csv(res, './Results/proportions/uncorrected/woProl_BisqueRNA.csv')
# 2.2.2. With correction before
res = Estimate.Proportions(CRC_bulk, list(data=references$corrected$data[, notprol_Tcells],
                                          metadata=references$corrected$metadata[notprol_Tcells,]), 'Deconv_cellTypes', 'BisqueRNA',
                           gene.markers=markers, use.overlap=FALSE, old.cpm=TRUE)
invisible(gc())
write.csv(res, './Results/proportions/corrected_before/woProl_BisqueRNA.csv')
# 2.2.3. With correction after
res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                          metadata=references$uncorrected$metadata[notprol_Tcells,]), 'Deconv_cellTypes', 'BisqueRNA',
                           gene.markers=markers, use.overlap=FALSE, old.cpm=TRUE)
invisible(gc())
res = correct_fractions_mRNABias(res, rna_correction_vec)
write.csv(res, './Results/proportions/corrected_after/woProl_BisqueRNA.csv')



# 3. BSeqSC
#    The method is run online, only missing getting results for the correction after, with and without
#    proliferative Tcells. All files need to be changed, as cell-type files should change from '_' to ' '.
# 3.1. With proliferative Tcells:
# 3.1.1. Uncorrected
res = read.csv('./Results/proportions/uncorrected/BSeqSC.csv', row.names=1)
rownames(res) = gsub('[_]', ' ', rownames(res))
write.csv(res, './Results/proportions/uncorrected/BSeqSC.csv')
# 3.1.2. Corected before
res = read.csv('./Results/proportions/corrected_before/BSeqSC.csv', row.names=1)
rownames(res) = gsub('[_]', ' ', rownames(res))
write.csv(res, './Results/proportions/corrected_before/BSeqSC.csv')
# 3.1.3. Corrected after
res = read.csv('./Results/proportions/uncorrected/BSeqSC.csv', row.names=1)
rownames(res) = gsub('[_]', ' ', rownames(res))
res = correct_fractions_mRNABias(res, rna_correction_vec)
write.csv(res, './Results/proportions/corrected_after/BSeqSC.csv')

# 3.2. Without proliferative Tcells:
# 3.2.1. Uncorrected
res = read.csv('./Results/proportions/uncorrected/woProl_BSeqSC.csv', row.names=1)
rownames(res) = gsub('[_]', ' ', rownames(res))
write.csv(res, './Results/proportions/uncorrected/woProl_BSeqSC.csv')
# 3.2.2. Corrected before
res = read.csv('./Results/proportions/corrected_before/woProl_BSeqSC.csv', row.names=1)
rownames(res) = gsub('[_]', ' ', rownames(res))
write.csv(res, './Results/proportions/corrected_before/woProl_BSeqSC.csv')
# 3.2.3. Corrected after
res = read.csv('./Results/proportions/uncorrected/woProl_BSeqSC.csv', row.names=1)
rownames(res) = gsub('[_]', ' ', rownames(res))
res = correct_fractions_mRNABias(res, rna_correction_vec)
write.csv(res, './Results/proportions/corrected_after/woProl_BSeqSC.csv')



# 4. CIBERSORTx
#    Default values: no batch correction; no permutations; 0.5 as minimum expression
#    The method is run online, only missing getting results for the correction after, with and without
#    proliferative Tcells
# 3.1. With proliferative Tcells:
# 3.1.1. Uncorrected
res = read.csv('./Results/proportions/uncorrected/CIBERSORTx.csv', row.names=1)
rownames(res) = gsub('[_]', ' ', rownames(res))
write.csv(res, './Results/proportions/uncorrected/CIBERSORTx.csv')
# 3.1.2. Corected before
res = read.csv('./Results/proportions/corrected_before/CIBERSORTx.csv', row.names=1)
rownames(res) = gsub('[_]', ' ', rownames(res))
write.csv(res, './Results/proportions/corrected_before/CIBERSORTx.csv')
# 3.1.3. Corrected after
res = read.csv('./Results/proportions/uncorrected/CIBERSORTx.csv', row.names=1)
rownames(res) = gsub('[_]', ' ', rownames(res))
res = correct_fractions_mRNABias(res, rna_correction_vec)
write.csv(res, './Results/proportions/corrected_after/CIBERSORTx.csv')

# 3.2. Without proliferative Tcells:
# 3.2.1. Uncorrected
res = read.csv('./Results/proportions/uncorrected/woProl_CIBERSORTx.csv', row.names=1)
rownames(res) = gsub('[_]', ' ', rownames(res))
write.csv(res, './Results/proportions/uncorrected/woProl_CIBERSORTx.csv')
# 3.2.2. Corrected before
res = read.csv('./Results/proportions/corrected_before/woProl_CIBERSORTx.csv', row.names=1)
rownames(res) = gsub('[_]', ' ', rownames(res))
write.csv(res, './Results/proportions/corrected_before/woProl_CIBERSORTx.csv')
# 3.2.3. Corrected after
res = read.csv('./Results/proportions/uncorrected/woProl_CIBERSORTx.csv', row.names=1)
rownames(res) = gsub('[_]', ' ', rownames(res))
res = correct_fractions_mRNABias(res, rna_correction_vec)
write.csv(res, './Results/proportions/corrected_after/woProl_CIBERSORTx.csv')



# 5. DigitalDLSorter

# 5.1. Run DigitalDLSorter, with proliferative cells
# 5.1.1. Without correction
ddls_res = Estimate.Proportions(CRC_bulk, references$uncorrected, 'Deconv_cellTypes', 'DigitalDLSorter', return.model=TRUE,
                                load_min.cells=10, simul_subset.cells=1100, simul_proportional=FALSE, simul_n.cells=1000,
                                simul_cell.types=c('NK cells', 'Proliferative Tcells'), bulk_num.bulk.samples=10000,
                                batch.size=128, threads=6)
invisible(gc())
write.csv(ddls_res$proportions, './Results/proportions/uncorrected/DigitalDLSorter.csv')
saveRDS(ddls_res$model_DDLS,
        '/home/scardoso/Documents/PhD/Tumour_Deconvolution/Results/DigitalDLSorter_models/uncorrected.Rdata')
# 5.1.2. With correction before
ddls_res = Estimate.Proportions(CRC_bulk, references$corrected, 'Deconv_cellTypes', 'DigitalDLSorter', return.model=TRUE,
                                load_min.cells=10, simul_subset.cells=1100, simul_proportional=FALSE, simul_n.cells=1000,
                                simul_cell.types=c('NK cells', 'Proliferative Tcells'), bulk_num.bulk.samples=10000,
                                batch.size=128, threads=6)
invisible(gc())
write.csv(ddls_res$proportions, './Results/proportions/corrected_before/DigitalDLSorter.csv')
saveRDS(ddls_res$model_DDLS,
        '/home/scardoso/Documents/PhD/Tumour_Deconvolution/Results/DigitalDLSorter_models/corrected.Rdata')
# 5.1.3. With correction after
res = read.csv('./Results/proportions/uncorrected/DigitalDLSorter.csv', row.names=1)
res = correct_fractions_mRNABias(res, rna_correction_vec)
invisible(gc())
write.csv(res, './Results/proportions/corrected_after/DigitalDLSorter.csv')

# 5.2. Run DigitalDLSorter, without proliferative cells
# 5.2.1. Without correction
ddls_res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                               metadata=references$uncorrected$metadata[notprol_Tcells,]),
                                'Deconv_cellTypes', 'DigitalDLSorter', return.model=TRUE,
                                load_min.cells=10, simul_subset.cells=1100, simul_proportional=FALSE, simul_n.cells=1000,
                                simul_cell.types=c('NK cells'), bulk_num.bulk.samples=10000,
                                batch.size=128, threads=6)
invisible(gc())
write.csv(ddls_res$proportions, './Results/proportions/uncorrected/woProl_DigitalDLSorter.csv')
saveRDS(ddls_res$model_DDLS,
        '/home/scardoso/Documents/PhD/Tumour_Deconvolution/Results/DigitalDLSorter_models/woProl_uncorrected.Rdata')
# 5.2.2. With correction before
ddls_res = Estimate.Proportions(CRC_bulk, list(data=references$corrected$data[, notprol_Tcells],
                                               metadata=references$corrected$metadata[notprol_Tcells,]),
                                'Deconv_cellTypes', 'DigitalDLSorter', return.model=TRUE,
                                load_min.cells=10, simul_subset.cells=1100, simul_proportional=FALSE, simul_n.cells=1000,
                                simul_cell.types=c('NK cells'), bulk_num.bulk.samples=10000,
                                batch.size=128, threads=6)
invisible(gc())
write.csv(ddls_res$proportions, './Results/proportions/corrected_before/woProl_DigitalDLSorter.csv')
saveRDS(ddls_res$model_DDLS,
        '/home/scardoso/Documents/PhD/Tumour_Deconvolution/Results/DigitalDLSorter_models/woProl_corrected.Rdata')
# 5.2.3. With correction after
res = read.csv('./Results/proportions/uncorrected/woProl_DigitalDLSorter.csv', row.names=1)
res = correct_fractions_mRNABias(res, rna_correction_vec)
invisible(gc())
write.csv(res, './Results/proportions/corrected_after/woProl_DigitalDLSorter.csv')



# 6. DWLS

# 6.1. Run DWLS, with proliferative Tcells
# 6.1.1. Without correction
res = Estimate.Proportions(CRC_bulk, references$uncorrected, 'Deconv_cellTypes', 'DWLS', 0.5, 0.01,
                           rna.bias=FALSE, bias.vec=NULL)
invisible(gc())
write.csv(res, './Results/proportions/uncorrected/DWLS.csv')
# 6.1.2. With correction after
res = Estimate.Proportions(CRC_bulk, references$uncorrected, 'Deconv_cellTypes', 'DWLS', 0.5, 0.01,
                           rna.bias=FALSE, bias.vec=NULL)
invisible(gc())
res = correct_fractions_mRNABias(res, rna_correction_vec)
write.csv(res, './Results/proportions/corrected_after/DWLS.csv')


# 6.2. Run DWLS, with proliferative Tcells
# 6.2.1. Without correction
res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                          metadata=references$uncorrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'DWLS', 0.5, 0.01, rna.bias=FALSE, bias.vec=NULL)
invisible(gc())
write.csv(res, './Results/proportions/uncorrected/woProl_DWLS.csv')
# 6.2.2. With correction after
res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                          metadata=references$uncorrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'DWLS', 0.5, 0.01, rna.bias=FALSE, bias.vec=NULL)
invisible(gc())
res = correct_fractions_mRNABias(res, rna_correction_vec)
write.csv(res, './Results/proportions/corrected_after/woProl_DWLS.csv')



# 7. MOMF

# 7.1. Run MOMF, with proliferative Tcells
# 7.1.1. Without correction
res = Estimate.Proportions(CRC_bulk, references$uncorrected, 'Deconv_cellTypes', 'MOMF', 'KL', 2,
                           1000, FALSE, NULL)
invisible(gc())
write.csv(res, './Results/proportions/uncorrected/MOMF.csv')
# 7.1.2. With before
res = Estimate.Proportions(CRC_bulk, references$corrected, 'Deconv_cellTypes', 'MOMF', 'KL', 2,
                           1000, FALSE, NULL)
invisible(gc())
write.csv(res, './Results/proportions/corrected_before/MOMF.csv')
# 7.1.3. With after
res = Estimate.Proportions(CRC_bulk, references$corrected, 'Deconv_cellTypes', 'MOMF', 'KL', 2,
                           1000, FALSE, NULL)
invisible(gc())
res = correct_fractions_mRNABias(res, rna_correction_vec)
write.csv(res, './Results/proportions/corrected_after/MOMF.csv')

# 7.2. Run MOMF, without proliferative Tcells
# 7.2.1. Without correction
res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                          metadata=references$uncorrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'MOMF', 'KL', 2, 1000, FALSE, NULL)
invisible(gc())
write.csv(res, './Results/proportions/uncorrected/woProl_MOMF.csv')
# 7.2.2. With before
res = Estimate.Proportions(CRC_bulk, list(data=references$corrected$data[, notprol_Tcells],
                                          metadata=references$corrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'MOMF', 'KL', 2, 1000, FALSE, NULL)
invisible(gc())
write.csv(res, './Results/proportions/corrected_before/woProl_MOMF.csv')
# 7.2.3. With after
res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                          metadata=references$uncorrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'MOMF', 'KL', 2, 1000, FALSE, NULL)
invisible(gc())
res = correct_fractions_mRNABias(res, rna_correction_vec)
write.csv(res, './Results/proportions/corrected_after/woProl_MOMF.csv')



# 8. MuSiC_woGrouping

# 8.1. Run MuSiC_woGrouping, with proliferative Tcells
res = Estimate.Proportions(CRC_bulk, references$uncorrected, 'Deconv_cellTypes', 'MuSiC',
                           music.method='without.Grouping', n.iterations=1000, center=FALSE,
                           normalize=FALSE, wo.method='Est.prop.weighted', wo.gene.markers=markers,
                           wo.bias.vec=rna_correction_vec)
invisible(gc())
write.csv(res, './Results/proportions/MuSiC_woGrouping.csv')

# 8.2. Run MuSiC_woGrouping, without proliferative Tcells
res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                          metadata=references$uncorrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'MuSiC', music.method='without.Grouping',
                           n.iterations=1000, center=FALSE, normalize=FALSE,
                           wo.method='Est.prop.weighted', wo.gene.markers=markers,
                           wo.bias.vec=rna_correction_vec)
invisible(gc())
write.csv(res, './Results/proportions/woProl_MuSiC_woGrouping.csv')



# 9. MuSiC_wGrouping

# 9.0. Get the clusters:
scReference = prepare_expressionSet_scReferences(references$uncorrected, 'Deconv_cellTypes')
design_matrix = MuSiC::music_basis(scReference, non.zero=TRUE, markers=markers,
                                   clusters='cellType', samples='sampleID', cell_size=NULL)
d = dist(t(log(design_matrix$Disgn.mtx + 1e-6)), method='euclidean')
hc = hclust(d, method='complete')
plot(hc, cex=0.6, hang=-1, main= 'Cluster log(Design Matrix)')
clusters = list(C1=c("Cancer cells", "Stromal cells", "Macro/mono Lineage"),
                C3=c("Bcells", "Regulatory CD4 Tcells", "CD4 Tcells", "Proliferative Tcells", "NK cells", "CD8 Tcells"))
clusters_woProl = list(C1=c("Cancer cells", "Stromal cells", "Macro/mono Lineage"),
                C3=c("Bcells", "Regulatory CD4 Tcells", "CD4 Tcells",  "NK cells", "CD8 Tcells"))

# 9.1. Run MuSiC_wGrouping, with proliferative cells
# 9.1.1. Without correction
res = Estimate.Proportions(CRC_bulk, references$uncorrected, 'Deconv_cellTypes', 'MuSiC', music.method='with.Grouping', n.iterations=1000,
                           center=FALSE, normalize=FALSE, w.celltypes.Cluster.List=clusters, w.markers.log2FC=.8, w.markers.pval=.01,
                           w.markers.only.Pos=TRUE)
invisible(gc())
write.csv(res, './Results/proportions/uncorrected/MuSiC_wGrouping.csv')
# 9.1.2. With before
res = Estimate.Proportions(CRC_bulk, references$corrected, 'Deconv_cellTypes', 'MuSiC', music.method='with.Grouping', n.iterations=1000,
                           center=FALSE, normalize=FALSE, w.celltypes.Cluster.List=clusters, w.markers.log2FC=.8, w.markers.pval=.01,
                           w.markers.only.Pos=TRUE)
invisible(gc())
write.csv(res, './Results/proportions/corrected_before/MuSiC_wGrouping.csv')
# 9.1.3. With after
res = Estimate.Proportions(CRC_bulk, references$uncorrected, 'Deconv_cellTypes', 'MuSiC', music.method='with.Grouping', n.iterations=1000,
                           center=FALSE, normalize=FALSE, w.celltypes.Cluster.List=clusters, w.markers.log2FC=.8, w.markers.pval=.01,
                           w.markers.only.Pos=TRUE)
invisible(gc())
res = correct_fractions_mRNABias(res, rna_correction_vec)
write.csv(res, './Results/proportions/corrected_after/MuSiC_wGrouping.csv')

# 9.2. Run MuSiC_wGrouping, with proliferative cells
# 9.2.1. Without correction
res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                          metadata=references$uncorrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'MuSiC', music.method='with.Grouping', n.iterations=1000, center=FALSE, normalize=FALSE,
                           w.celltypes.Cluster.List=clusters_woProl, w.markers.log2FC=.8, w.markers.pval=.01, w.markers.only.Pos=TRUE)
invisible(gc())
write.csv(res, './Results/proportions/uncorrected/woProl_MuSiC_wGrouping.csv')
# 9.2.2. With before
res = Estimate.Proportions(CRC_bulk, list(data=references$corrected$data[, notprol_Tcells],
                                          metadata=references$corrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'MuSiC', music.method='with.Grouping', n.iterations=1000, center=FALSE, normalize=FALSE, 
                          w.celltypes.Cluster.List=clusters_woProl, w.markers.log2FC=.8, w.markers.pval=.01, w.markers.only.Pos=TRUE)
invisible(gc())
write.csv(res, './Results/proportions/corrected_before/woProl_MuSiC_wGrouping.csv')
# 9.2.3. With after
res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                          metadata=references$uncorrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'MuSiC', music.method='with.Grouping', n.iterations=1000, center=FALSE, normalize=FALSE,
                           w.celltypes.Cluster.List=clusters_woProl, w.markers.log2FC=.8, w.markers.pval=.01, w.markers.only.Pos=TRUE)
invisible(gc())
res = correct_fractions_mRNABias(res, rna_correction_vec)
write.csv(res, './Results/proportions/corrected_after/woProl_MuSiC_wGrouping.csv')



# 10. Scaden

# 10.0. Scaden changes non-alphanumeric caracters into spaces:
cell_types_scaden = cell_types
names(cell_types_scaden) = cell_types
cell_types_scaden[3] = 'Macro mono Lineage'

# 10.1. Run Scaden, with proliferative cells
# 10.1.1. Without correction
res = Estimate.Proportions(CRC_bulk_file_scaden, scaden_train_files$uncorrected$wProl, 'Deconv_cellTypes',
                           'Scaden', dir.results=tempdir(), datasets='', min.expression=0.1,
                           batch.size=128L, learning.rate=0.0001, n.steps=1000L)
rownames(res)[3] = 'Macro/mono Lineage'
invisible(gc())
write.csv(res, './Results/proportions/uncorrected/Scaden.csv')
# 10.1.2. With correction before
res = Estimate.Proportions(CRC_bulk_file_scaden, scaden_train_files$corrected$wProl, 'Deconv_cellTypes',
                           'Scaden', dir.results=tempdir(), datasets='', min.expression=0.1,
                           batch.size=128L, learning.rate=0.0001, n.steps=1000L)
rownames(res)[3] = 'Macro/mono Lineage'
invisible(gc())
write.csv(res, './Results/proportions/corrected_before/Scaden.csv')
# 10.1.3. With correction after
res = Estimate.Proportions(CRC_bulk_file_scaden, scaden_train_files$uncorrected$wProl, 'Deconv_cellTypes',
                           'Scaden', dir.results=tempdir(), datasets='', min.expression=0.1,
                           batch.size=128L, learning.rate=0.0001, n.steps=1000L)
rownames(res)[3] = 'Macro/mono Lineage'
res = correct_fractions_mRNABias(res, rna_correction_vec)
invisible(gc())
write.csv(res, './Results/proportions/corrected_after/Scaden.csv')

# 10.2. Run Scaden, without proliferative cells
# 10.2.1. Without correction
res = Estimate.Proportions(CRC_bulk_file_scaden, scaden_train_files$uncorrected$woProl, 'Deconv_cellTypes',
                           'Scaden', dir.results=tempdir(), datasets='', min.expression=0.1,
                           batch.size=128L, learning.rate=0.0001, n.steps=1000L)
rownames(res)[3] = 'Macro/mono Lineage'
invisible(gc())
write.csv(res, './Results/proportions/uncorrected/woProl_Scaden.csv')
# 10.2.2. With correction before
res = Estimate.Proportions(CRC_bulk_file_scaden, scaden_train_files$corrected$woProl, 'Deconv_cellTypes',
                           'Scaden', dir.results=tempdir(), datasets='', min.expression=0.1,
                           batch.size=128L, learning.rate=0.0001, n.steps=1000L)
rownames(res)[3] = 'Macro/mono Lineage'
invisible(gc())
write.csv(res, './Results/proportions/corrected_before/woProl_Scaden.csv')
# 10.2.3. With correction after
res = Estimate.Proportions(CRC_bulk_file_scaden, scaden_train_files$uncorrected$woProl, 'Deconv_cellTypes',
                           'Scaden', dir.results=tempdir(), datasets='', min.expression=0.1,
                           batch.size=128L, learning.rate=0.0001, n.steps=1000L)
rownames(res)[3] = 'Macro/mono Lineage'
res = correct_fractions_mRNABias(res, rna_correction_vec)
invisible(gc())
write.csv(res, './Results/proportions/corrected_after/woProl_Scaden.csv')


# 11. SCDC

# 11.1. Run SCDC, with proliferative Tcells
# 11.1.1. Without correction
res = Estimate.Proportions(CRC_bulk, references$uncorrected, 'Deconv_cellTypes', 'SCDC',
                           scdc.method='multiple', n.iterations=1000, m.use.gridSearch=FALSE,
                           m.search.length=0.05, m.dataset.var='Dataset', m.weights.method='LAD')
invisible(gc())
write.csv(res, './Results/proportions/uncorrected/SCDC.csv')
# 11.1.2. With correction before
res = Estimate.Proportions(CRC_bulk, references$corrected, 'Deconv_cellTypes', 'SCDC',
                           scdc.method='multiple', n.iterations=1000, m.use.gridSearch=FALSE,
                           m.search.length=0.05, m.dataset.var='Dataset', m.weights.method='LAD')
invisible(gc())
write.csv(res, './Results/proportions/corrected_before/SCDC.csv')
# 11.1.2. With correction after
res = Estimate.Proportions(CRC_bulk, references$uncorrected, 'Deconv_cellTypes', 'SCDC',
                           scdc.method='multiple', n.iterations=1000, m.use.gridSearch=FALSE,
                           m.search.length=0.05, m.dataset.var='Dataset', m.weights.method='LAD')
invisible(gc())
res = correct_fractions_mRNABias(res, rna_correction_vec)
write.csv(res, './Results/proportions/corrected_after/SCDC.csv')

# 11.2. Run SCDC, without proliferative Tcells
# 11.2.1. Without correction
res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                          metadata=references$uncorrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'SCDC', scdc.method='multiple', n.iterations=1000, m.use.gridSearch=FALSE,
                           m.search.length=0.05, m.dataset.var='Dataset', m.weights.method='LAD')
invisible(gc())
write.csv(res, './Results/proportions/uncorrected/woProl_SCDC.csv')
# 11.2.2. With correction before
res = Estimate.Proportions(CRC_bulk, list(data=references$corrected$data[, notprol_Tcells],
                                          metadata=references$corrected$metadata[notprol_Tcells,]), 'Deconv_cellTypes', 'SCDC',
                           scdc.method='multiple', n.iterations=1000, m.use.gridSearch=FALSE,
                           m.search.length=0.05, m.dataset.var='Dataset', m.weights.method='LAD')
invisible(gc())
write.csv(res, './Results/proportions/corrected_before/woProl_SCDC.csv')
# 11.2.2. With correction after
res = Estimate.Proportions(CRC_bulk, list(data=references$uncorrected$data[, notprol_Tcells],
                                          metadata=references$uncorrected$metadata[notprol_Tcells,]),
                           'Deconv_cellTypes', 'SCDC', scdc.method='multiple', n.iterations=1000, m.use.gridSearch=FALSE,
                           m.search.length=0.05, m.dataset.var='Dataset', m.weights.method='LAD')
invisible(gc())
res = correct_fractions_mRNABias(res, rna_correction_vec)
write.csv(res, './Results/proportions/corrected_after/woProl_SCDC.csv')




