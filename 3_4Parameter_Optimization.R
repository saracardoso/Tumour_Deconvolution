code_dir = './code/R'
for(file in list.files(code_dir, full.names=TRUE)) source(file)

# methods = c('AutoGeneS', 'BisqueRNA', 'BSeqSC', 'CIBERSORTx', 'DigitalDLSorter', 'DWLS', 'MOMF', 'MuSiC_woGrouping', 'MuSiC_wGrouping',
#             'Scaden', 'SCDC')
methods = c('AutoGeneS', 'BisqueRNA', 'BSeqSC', 'CIBERSORTx', 'DWLS', 'MOMF', 'MuSiC_woGrouping', 'MuSiC_wGrouping', 'Scaden', 'SCDC')
cell_types = c('Cancer cells', 'Stromal cells', 'Anti-Inflammatory macro/mono', 'Pro-Inflammatory macro/mono', 'Bcells', 'CD4 Tcells',
               'Regulatory CD4 Tcells', 'CD8 Tcells', 'Proliferative Tcells', 'NK cells', 'Other cells')



# ---
# - Read bulk and reference files for test optimization
# ---

# 1. Read bulk data:
test = readRDS('./Data/bulk/CRC/NICs/train_test_sets/test.Rdata')
test_file_scaden = './Data/bulk/CRC/NICs/train_test_sets/test.txt'



# 2. Read reference files:
# 2.1. Reference counts:
references = list()
references$uncorrected = read_scRNAseq_reference('./Data/CRC_reference/UB_matrix_CSV.csv', './Data/CRC_reference/metadata.csv')
invisible(gc())
references$corrected = read_scRNAseq_reference('./Data/CRC_reference/B_matrix_CSV.csv', './Data/CRC_reference/metadata.csv')
invisible(gc())

# 2.2. Gene Markers:
markers =  jsonlite::read_json('./Data/CRC_reference/markers.json', simplifyVector=T)

# 2.3. Scaden train file:
scaden_train_files = list()
scaden_train_files$corrected = './Data/CRC_reference/UB_Scaden.h5ad'
scaden_train_files$uncorrected = './Data/CRC_reference/B_Scaden.h5ad'

# 2.4. RNA bias correction vector
rna_correction_df = read.csv('./Data/CRC_reference/bias_values.csv', row.names=1)
rna_correction_vec = rna_correction_df[,1]
names(rna_correction_vec) = rownames(rna_correction_df)





# ---
# - Read methods' parameters tested in train optimization
# ---
methods_parameters = list()
for(method in methods){
  methods_parameters[[method]] = read.csv(paste('./Results/1_Parameter_Optimization/parameters/', method, '.csv', sep=''),
                                          row.names=1, stringsAsFactors=FALSE)
}





# ---
# - Get best combinations from each method, calculated using the training set
# ---
best_combinations = c()
for(method in methods){
  if(method == 'DigitalDLSorter') next
  train_results = read.csv(paste('./Results/1_Parameter_Optimization/train/', method, '.csv', sep=''), header=TRUE, row.names=1)
  best_combinations = c(best_combinations, names(which.min(train_results['MEAN',])))
}
names(best_combinations) = methods





# ---
# - Parameter optimization: Test set
# ---



# 1. AutoGeneS

# 1.1. Run test
best_combo = best_combinations['AutoGeneS']
combo = methods_parameters$AutoGeneS[best_combo,]
gene_markers = NULL
if(combo$gene.markers=='YES') gene_markers = markers
w = list(-1,1)
if(combo$weights=='c'){
  w = list(-1,0)
} else if(combo$weights=='d'){
  w = list(0,-1)
}
if(combo$n.genes == 'all'){
  if(!is.null(gene_markers)) n_genes = as.integer(length(unique(unlist(gene_markers)))) # all gene markers
  else n_genes = as.integer(800)
} else{
  n_genes = as.integer(combo$n.genes)
}
# Estimated proportions:
res_test = Estimate.Proportions(test, references$uncorrected, 'Deconv_cellTypes', 'AutoGeneS', model=combo$model, gene.markers=gene_markers,
                                mode=combo$mode, n.iterations=as.integer(combo$n.iterations), n.genes=as.integer(combo$n.genes), weights=w)
res_test[res_test<0] = 0
invisible(gc())

# 1.2. Save estimated proportions:
write.csv(res_test, './Results/1_Parameter_Optimization/test/proportions/AutoGeneS.csv')



# 2. BisqueRNA

# 2.1. Run test
best_combo = best_combinations['BisqueRNA']
combo = methods_parameters$BisqueRNA[tail(strsplit(best_combo, '_')[[1]], n=1),]
gene_markers = NULL
if(combo$gene.markers=='YES') gene_markers = markers
correction = 'uncorrected'
if(strsplit(best_combo, '_')[[1]][1] == 'corrected.Before') correction = 'corrected'
res_test = Estimate.Proportions(test, references[[correction]], 'Deconv_cellTypes', 'BisqueRNA', gene.markers=gene_markers,
                                use.overlap=FALSE, old.cpm=combo$old.cpm)
if(strsplit(best_combo, '_')[[1]][1] == 'corrected.After') res_test = correct_fractions_mRNABias(res_test, rna_correction_vec)
invisible(gc())

# 2.2. Save estimated proportions:
write.csv(res_test, './Results/1_Parameter_Optimization/test/proportions/BisqueRNA.csv')



# 3. BSeqSC
best_combo = best_combinations['BSeqSC']
best_combo
combo = methods_parameters$BSeqSC[tail(strsplit(best_combo, '_')[[1]], n=1),]
combo

# 3.1. Run test
res_test = t(read.csv('./Results/1_Parameter_Optimization/test/BSeqSC_web/BSeqSC_test.csv', row.names=1))
res_test = res_test[1:11,]
rownames(res_test) = gsub('[_]', ' ', rownames(res_test))
rownames(res_test)[3:4] = c('Anti-Inflammatory macro/mono', 'Pro-Inflammatory macro/mono')

# 3.2. Save estimated of proportions:
write.csv(res_test, './Results/1_Parameter_Optimization/test/proportions/BSeqSC.csv')



# 4. CIBERSORTx
best_combo = best_combinations['CIBERSORTx']
best_combo
combo = methods_parameters$CIBERSORTx[tail(strsplit(best_combo, '_')[[1]], n=1),]
combo

# 4.1. Run test
res_test = t(read.csv('./Results/1_Parameter_Optimization/test/CIBERSORTx_web/CIBERSORTx_test.csv', row.names=1))
res_test = res_test[1:11,]
rownames(res_test) = gsub('[_]', ' ', rownames(res_test))
rownames(res_test)[3:4] = c('Anti-Inflammatory macro/mono', 'Pro-Inflammatory macro/mono')

# 4.2. Save estimated of proportions:
write.csv(res_test, './Results/1_Parameter_Optimization/test/proportions/CIBERSORTx.csv')



# 5. DigitalDLSorter

# 5.1. Run test
best_combo = best_combinations['DigitalDLSorter']
combo = methods_parameters$DWLS[tail(strsplit(best_combo, '_')[[1]], n=1),]
DDLSorter_model = read.RDS('./Results/1_Parameter_Optimization/train/DigitalDLSorter_utils/trained_model.Rdata')

res_test = Estimate.Proportions(test, DDLSorter_model, 'Deconv_cellTypes', 'DigitalDLSorter', pipeline='predict',
                                predict_normalize=combo$predict_normalize, batch.size=128, threads=6)
invisible(gc())

# 5.2. Save estimated of proportions:
write.csv(res_test$proportions, './Results/1_Parameter_Optimization/test/proportions/DigitalDLSorter.csv')



# 6. DWLS

# 6.1. Run test
best_combo = best_combinations['DWLS']
combo = methods_parameters$DWLS[tail(strsplit(best_combo, '_')[[1]], n=1),]
correction = 'uncorrected'
if(strsplit(best_combo, '_')[[1]][1] == 'corrected.Before') correction = 'corrected'
res_test = Estimate.Proportions(test, references[[correction]], 'Deconv_cellTypes', 'DWLS', combo$diff.cutoff, combo$pval.cutoff,
                                rna.bias=FALSE, bias.vec=NULL)
if(strsplit(best_combo, '_')[[1]][1] == 'corrected.After') res_test = correct_fractions_mRNABias(res_test, rna_correction_vec)
res_test[res_test<0] = 0
invisible(gc())

# 6.2. Save estimated proportions:
write.csv(res_test, './Results/1_Parameter_Optimization/test/proportions/DWLS.csv')



# 7. MOMF

# 7.1. Run test
best_combo = best_combinations['MOMF']
combo = methods_parameters$MOMF[tail(strsplit(best_combo, '_')[[1]], n=1),]
correction = 'uncorrected'
if(strsplit(best_combo, '_')[[1]][1] == 'corrected.Before') correction = 'corrected'
res_test = Estimate.Proportions(test, references[[correction]], 'Deconv_cellTypes', 'MOMF', combo$method, combo$rho, combo$num.iter,
                                FALSE, NULL)
if(strsplit(best_combo, '_')[[1]][1] == 'corrected.After') res_test = correct_fractions_mRNABias(res_test, rna_correction_vec)
invisible(gc())

# 7.2. Save estimated proportions:
write.csv(res_test, './Results/1_Parameter_Optimization/test/proportions/MOMF.csv')



# 8. MuSiC_woGrouping

# 8.1. Run test
best_combo = best_combinations['MuSiC_woGrouping']
combo = methods_parameters$MuSiC_woGrouping[best_combo,]
gene_markers = NULL
if(combo$wo.gene.markers=='YES') gene_markers = markers
method = 'Est.prop.weighted'
if(combo$wo.method=='all genes') method = 'Est.prop.allgene'
res_test = Estimate.Proportions(test, references$uncorrected, 'Deconv_cellTypes', 'MuSiC', music.method='without.Grouping', n.iterations=combo$n.iterations,
                                center=combo$center, normalize=combo$normalise, wo.method=method,
                                wo.gene.markers=gene_markers, wo.bias.vec=rna_correction_vec)
invisible(gc())

# 8.2. Save estimated proportions:
write.csv(res_test, './Results/1_Parameter_Optimization/test/proportions/MuSiC_woGrouping.csv')



# 9. MuSiC_wGrouping

# 9.1. Run test
best_combo = best_combinations['MuSiC_wGrouping']
combo = methods_parameters$MuSiC_wGrouping[tail(strsplit(best_combo, '_')[[1]], n=1),]
clusters = list(C1=c("Cancer cells", "Stromal cells"), C2=c("Anti-Inflammatory macro/mono", "Pro-Inflammatory macro/mono", "Other cells"),
                C3=c("Bcells", "Regulatory CD4 Tcells", "CD4 Tcells", "Proliferative Tcells", "NK cells", "CD8 Tcells"))
correction = 'uncorrected'
if(strsplit(best_combo, '_')[[1]][1] == 'corrected.Before') correction = 'corrected'
res_test = Estimate.Proportions(test, references[[correction]], 'Deconv_cellTypes', 'MuSiC', music.method='with.Grouping',
                                n.iterations=combo$n.iterations, center=combo$center, normalize=combo$normalise,
                                w.celltypes.Cluster.List=clusters, w.markers.log2FC=.8,
                                w.markers.pval=.01, w.markers.only.Pos=TRUE)
if(strsplit(best_combo, '_')[[1]][1] == 'corrected.After') res_test = correct_fractions_mRNABias(res_test, rna_correction_vec)
invisible(gc())

# 9.2. Save estimated proportions:
write.csv(res_test, './Results/1_Parameter_Optimization/test/proportions/MuSiC_wGrouping.csv')



# 10. Scaden
cell_types_scaden = cell_types
names(cell_types_scaden) = cell_types
cell_types_scaden[3:4] = c('Anti Inflammatory macro mono', 'Pro Inflammatory macro mono')
rna_correction_vec_scaden = rna_correction_vec
names(rna_correction_vec_scaden) = cell_types_scaden[match(names(rna_correction_vec), names(cell_types_scaden))]

# 10.1. Run test
best_combo = best_combinations['Scaden']
combo = methods_parameters$Scaden[tail(strsplit(best_combo, '_')[[1]], n=1),]
correction = 'uncorrected'
if(strsplit(best_combo, '_')[[1]][1] == 'corrected.Before') correction = 'corrected'
res_test = Estimate.Proportions(test_file_scaden, scaden_train_files[[correction]], 'Deconv_cellTypes', 'Scaden', dir.results=tempdir(), datasets='',
                                min.expression=combo$min.expression, batch.size=as.integer(combo$batch.size),
                                learning.rate=combo$learning.rate, n.steps=as.integer(combo$n.steps))
if(strsplit(best_combo, '_')[[1]][1] == 'corrected.After') res_test = correct_fractions_mRNABias(res_test, rna_correction_vec_scaden)
rownames(res_test) = names(cell_types_scaden)
invisible(gc())

# 10.2. Save estimated proportions:
write.csv(res_test, './Results/1_Parameter_Optimization/test/proportions/Scaden.csv')



# 11. SCDC
best_combo = best_combinations['SCDC']
combo = methods_parameters$SCDC[tail(strsplit(best_combo, '_')[[1]], n=1),]
weights_method = strsplit(best_combo, '[.]')[[1]][1]
correction = 'uncorrected'
if(strsplit(best_combo, '_')[[1]][1] == 'corrected.Before') correction = 'corrected'
res_test = Estimate.Proportions(test, references[[correction]], 'Deconv_cellTypes', 'SCDC', scdc.method='multiple',
                                n.iterations=combo$m.n.iterations, m.use.gridSearch=FALSE, m.search.length=combo$m.search.length,
                                m.dataset.var='Dataset', m.weights.method=weights_method)[[weights_method]]
if(strsplit(best_combo, '_')[[1]][1] == paste(weights_method, 'corrected.After', sep='.'))
  res_test = correct_fractions_mRNABias(res_test, rna_correction_vec)
invisible(gc())

# 11.2. Save estimated proportions:
write.csv(res_test, './Results/1_Parameter_Optimization/test/proportions/SCDC.csv')


