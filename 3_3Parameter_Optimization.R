code_dir = './code/R'
for(file in list.files(code_dir, full.names=TRUE)) source(file)

methods = c('AutoGeneS', 'BisqueRNA', 'BSeqSC', 'CIBERSORTx', 'DigitalDLSorter', 'DWLS', 'MOMF', 'MuSiC_woGrouping', 'MuSiC_wGrouping', 'Scaden', 'SCDC')
cell_types = c('Cancer cells', 'Stromal cells', 'Anti-Inflammatory macro/mono', 'Pro-Inflammatory macro/mono', 'Bcells', 'CD4 Tcells',
               'Regulatory CD4 Tcells', 'CD8 Tcells', 'Proliferative Tcells', 'NK cells', 'Other cells')





# ---
# - Read bulk and reference files for optimization
# ---

# 1. Read bulk data:
train = readRDS('./Data/bulk/CRC/NICs/train_test_sets/train.Rdata')
train_bulk_file_scaden = './Data/bulk/CRC/NICs/train_test_sets/train.txt'


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
scaden_train_files$corrected = './Data/CRC_reference/B_Scaden.h5ad'
scaden_train_files$uncorrected = './Data/CRC_reference/UB_Scaden.h5ad'

# 2.4. RNA bias correction vector
rna_correction_df = read.csv('./Data/CRC_reference/bias_values.csv', row.names=1)
rna_correction_vec = rna_correction_df[,1]
names(rna_correction_vec) = rownames(rna_correction_df)





# ---
# - Read methods' parameters for optimization
# ---
methods_parameters = list()
for(method in methods){
  methods_parameters[[method]] = read.csv(paste('./Results/1_Parameter_Optimization/parameters/', method, '.csv', sep=''),
                                          row.names=1, stringsAsFactors=FALSE)
}





# ---
# - Parameter optimization: Train set
# ---

# 1. AutoGeneS
train_rmses = c()
train_rmses_samp = c()
options = c()

# 1.1. Training:
n_params = dim(methods_parameters$AutoGeneS)[1]
pb = txtProgressBar(min=0, max=n_params, style=3, width=n_params, char='=')
init = numeric(n_params)
end = numeric(n_params)
for(combo_i in 1:n_params){
  init[combo_i] <- Sys.time()
  combo_name = rownames(methods_parameters$AutoGeneS)[combo_i]
  options = c(options, combo_name)
  message(paste(combo_name, '|', combo_i, '/', length(init), sep=' '))
  combo = methods_parameters$AutoGeneS[combo_name, ]
  gene_markers = NULL
  if(combo$gene.markers=='YES') gene_markers = markers
  w = list(-1,1)
  if(combo$weights=='c') w = list(-1,0)
  else if(combo$weights=='d') w = list(0,-1)
  if(combo$n.genes == 'all')
    if(!is.null(gene_markers)) n_genes = as.integer(length(unique(unlist(gene_markers)))) # all gene markers
    else n_genes = as.integer(800)
  else n_genes = as.integer(combo$n.genes)
  res = Estimate.Proportions(train, references$uncorrected, 'Deconv_cellTypes', 'AutoGeneS', model=combo$model, gene.markers=gene_markers,
                             mode=combo$mode, n.iterations=as.integer(combo$n.iterations), n.genes=n_genes, weights=w)
  invisible(gc())
  if(correction=='corrected.After') res = correct_fractions_mRNABias(res, rna_correction_vec)
  rmses_samp = c()
  for(samp in colnames(res)) rmses_samp = c(rmses_samp, Metrics::rmse(unlist(train$metadata[cell_types, samp]),
                                                                      unlist(res[cell_types, samp])))
  train_rmses_samp = cbind(train_rmses_samp, rmses_samp)
  train_rmses = c(train_rmses, mean(rmses_samp))
  message('')
  end[ti] = Sys.time()
  setTxtProgressBar(pb, ti)
  time = round(lubridate::seconds_to_period(sum(end - init)), 0)
  est = n_params * (mean(end[end != 0] - init[init != 0])) - sum(end - init)
  remainining <- round(lubridate::seconds_to_period(est), 0)
  cat(paste(" // Execution time:", time, " // Estimated time remaining:", remainining), "\n")
  invisible(gc())
}
names(train_rmses) = options
colnames(train_rmses_samp) = options
rownames(train_rmses_samp) = colnames(res)

# 1.2. Save train results in csv:
full_train_results = rbind(train_rmses_samp, train_rmses)
rownames(full_train_results)[dim(full_train_results)[1]] = 'MEAN'
write.csv(full_train_results, './Results/1_Parameter_Optimization/train/AutoGeneS.csv')



# 2. BisqueRNA
train_rmses = c()
train_rmses_samp = c()
options = c()

# 2.1. Training:
n_params = dim(methods_parameters$BisqueRNA)[1]
pb = txtProgressBar(min=0, max=n_params * 3, style=3, width=n_params * 3, char='=')
init = numeric(n_params * 3)
end = numeric(n_params * 3)
cr = 1
for(correction in c('corrected.Before', 'corrected.After', 'uncorrected')){
  if(correction=='corrected.Before') reference = references[['corrected']]
  else reference = references[['uncorrected']]
  for(combo_i in 1:n_params){
    if(cr==1) ti = combo_i
    if(cr==2) ti = n_params + combo_i
    if(cr==3) ti = (n_params*2) + combo_i
    init[ti] <- Sys.time()
    combo_name = rownames(methods_parameters$BisqueRNA)[combo_i]
    options = c(options, paste(correction, combo_name, sep='_'))
    message(paste(combo_name, '|', ti, '/', length(init), sep=' '))
    combo = methods_parameters$BisqueRNA[combo_name, ]
    gene_markers = NULL
    if(combo$gene.markers=='YES') gene_markers = markers
    res = Estimate.Proportions(train, reference, 'Deconv_cellTypes', 'BisqueRNA', gene.markers=gene_markers, use.overlap=FALSE, old.cpm=combo$old.cpm)
    if(correction=='corrected.After') res = correct_fractions_mRNABias(res, rna_correction_vec)
    rmses_samp = c()
    for(samp in colnames(res)) rmses_samp = c(rmses_samp, Metrics::rmse(unlist(train$metadata[cell_types, samp]),
                                                                        unlist(res[cell_types, samp])))
    train_rmses_samp = cbind(train_rmses_samp, rmses_samp)
    train_rmses = c(train_rmses, mean(rmses_samp))
    message('')
    end[ti] = Sys.time()
    setTxtProgressBar(pb, ti)
    time = round(lubridate::seconds_to_period(sum(end - init)), 0)
    est = n_params * 3 * (mean(end[end != 0] - init[init != 0])) - sum(end - init)
    remainining <- round(lubridate::seconds_to_period(est), 0)
    cat(paste(" // Execution time:", time, " // Estimated time remaining:", remainining), "\n")
    invisible(gc())
  }
  remove(reference)
  invisible(gc())
  cr = cr+1
}
names(train_rmses) = options
time = end - init
names(time) = options
colnames(train_rmses_samp) = options
rownames(train_rmses_samp) = colnames(res)

# 2.2. Save train results in csv:
full_train_results = rbind(train_rmses_samp, train_rmses, time)
nrows = dim(full_train_results)[1]
rownames(full_train_results)[(nrows-1):nrows] = c('MEAN', 'TIME')
write.csv(full_train_results, './Results/1_Parameter_Optimization/train/BisqueRNA.csv')



# 3. BSeqSC
# This method is run through the CIBERSORTx's website. As such, here we only have the organization of the results into the appropriate format
train_files_BSeqSC = list.files('./Results/1_Parameter_Optimization/train/BSeqSC_web')
train_rmses = c()
train_rmses_samp = c()
options = c()
for(res_file in train_files_BSeqSC){
  combo_name = gsub('[.]csv', '', res_file)
  res = t(read.csv(paste('./Results/1_Parameter_Optimization/train/BSeqSC_web', res_file, sep='/'), row.names=1, check.names=FALSE))
  rownames(res) = gsub('[_]', ' ', rownames(res))
  res = res[cell_types,]
  if(length(grep('uncorrected', combo_name))==1){
    for(correction in c('uncorrected', 'corrected.After')){
      if(correction=='corrected.After'){
        combo_name = gsub('uncorrected', 'corrected.After', combo_name)
        res = correct_fractions_mRNABias(res, rna_correction_vec)
      }
      rmses_samp = c()
      for(samp in colnames(res)) rmses_samp = c(rmses_samp, Metrics::rmse(unlist(train$metadata[cell_types, samp]),
                                                                          unlist(res[cell_types, samp])))
      options = c(options, combo_name)
      train_rmses_samp = cbind(train_rmses_samp, rmses_samp)
      train_rmses = c(train_rmses, mean(rmses_samp))
    }
  }
  else{
    rmses_samp = c()
    for(samp in colnames(res)) rmses_samp = c(rmses_samp, Metrics::rmse(unlist(train$metadata[cell_types, samp]),
                                                                        unlist(res[cell_types, samp])))
    options = c(options, combo_name)
    train_rmses_samp = cbind(train_rmses_samp, rmses_samp)
    train_rmses = c(train_rmses, mean(rmses_samp))
  }
}
names(train_rmses) = options
colnames(train_rmses_samp) = options
rownames(train_rmses_samp) = colnames(res)
# 3.1. Save train results in csv:
full_train_results = rbind(train_rmses_samp, train_rmses)
rownames(full_train_results)[dim(full_train_results)[1]] = 'MEAN'
write.csv(full_train_results, './Results/1_Parameter_Optimization/train/BSeqSC.csv')



# 4. CIBERSORTx
# This method is run through their website. As such, here we only have the organization of the results into the appropriated format
train_files_cibersortx = list.files('./Results/1_Parameter_Optimization/train/CIBERSORTx_web')
train_rmses = c()
train_rmses_samp = c()
options = c()
for(res_file in train_files_cibersortx){
  combo_name = gsub('[.]csv', '', res_file)
  res = t(read.csv(paste('./Results/1_Parameter_Optimization/train/CIBERSORTx_web',  res_file, sep='/'), row.names=1, check.names=FALSE))
  rownames(res) = gsub('[_]', ' ', rownames(res))
  res = res[cell_types,]
  if(length(grep('uncorrected', combo_name))==1){
    for(correction in c('uncorrected', 'corrected.After')){
      if(correction=='corrected.After'){
        combo_name = gsub('uncorrected', 'corrected.After', combo_name)
        res = correct_fractions_mRNABias(res, rna_correction_vec)
      } 
      rmses_samp = c()
      for(samp in colnames(res)) rmses_samp = c(rmses_samp, Metrics::rmse(unlist(train$metadata[cell_types, samp]),
                                                                          unlist(res[cell_types, samp])))
      options = c(options, combo_name)
      train_rmses_samp = cbind(train_rmses_samp, rmses_samp)
      train_rmses = c(train_rmses, mean(rmses_samp))
    }
  }
  else{
    rmses_samp = c()
    for(samp in colnames(res)) rmses_samp = c(rmses_samp, Metrics::rmse(unlist(train$metadata[cell_types, samp]),
                                                                        unlist(res[cell_types, samp])))
    options = c(options, combo_name)
    train_rmses_samp = cbind(train_rmses_samp, rmses_samp)
    train_rmses = c(train_rmses, mean(rmses_samp))
  }
}
names(train_rmses) = options
colnames(train_rmses_samp) = options
rownames(train_rmses_samp) = colnames(res)
# 4.1. Save train results in csv:
full_train_results = rbind(train_rmses_samp, train_rmses)
rownames(full_train_results)[dim(full_train_results)[1]] = 'MEAN'
write.csv(full_train_results, './Results/1_Parameter_Optimization/train/CIBERSORTx.csv')



# 5. DigitalDLSorter

# 5.1. Simulate cells:
DDLS_objects = list()
cell_simul_options = expand.grid(c(0, 1), c('limit_1000'), c(110, 1100), c('YES', 'NO'), stringsAsFactors=F)
rownames(cell_simul_options) = apply(expand.grid(c('A1', 'A2'), c('B2'), c('C1', 'C2'), c('D1', 'D2'),
                                                 stringsAsFactors=F), 1, paste, collapse='')
colnames(cell_simul_options) = c('load_min.counts', 'simulation', 'simul_subset.cells', 'simul_proportional')
n_params_sim = dim(cell_simul_options)[1]
pb_sim = txtProgressBar(min=0, max=n_params_sim * 2, style=3, width=200, char='=')
init_sim = numeric(n_params_sim * 2)
end_sim = numeric(n_params_sim * 2)
cr=1
for(ref in c('corrected', 'uncorrected')){
  DDLS_objects[[ref]] = list()
  for(combo_i in 1:n_params_sim){
    if(cr==1) ti = combo_i
    if(cr==2) ti = n_params + combo_i
    init_sim[ti] <- Sys.time()
    combo_name = rownames(cell_simul_options)[combo_i]
    message(paste(combo_name, '|', ti, '/', length(init_sim) * 2, sep=' '))
    combo = cell_simul_options[combo_name, ]
    
    simul_proportional = FALSE
    if(combo$simul_proportional=='YES') simul_proportional = TRUE
    DDLSorter_obj = prepare_singlecellexperiment_scReferences(references[[ref]], 'Deconv_cellTypes', 10, combo$load_min.counts, 'DDLS')
    DDLSorter_obj = simulate_cells(DDLSorter_obj, 1000, c('NK cells', 'Proliferative Tcells'), combo$simul_subset.cells,
                                   simul_proportional, 5, 'limited')
    saveRDS(DDLSorter_obj, paste('./Results/1_Parameter_Optimization/train/DigitalDLSorter_utils/simulate_cells/', ref, '_', combo_name,
                                 '.Rdata', sep=''))
    remove(DDLSorter_obj)
    invisible(gc())
    
    message('')
    end_sim[ti] = Sys.time()
    setTxtProgressBar(pb_sim, ti)
    time = round(lubridate::seconds_to_period(sum(end_sim - init_sim)), 0)
    est = n_params_sim * 2 * (mean(end_sim[end_sim != 0] - init_sim[init_sim != 0])) - sum(end_sim - init_sim)
    remainining <- round(lubridate::seconds_to_period(est), 0)
    cat(paste(" // Execution time:", time, " // Estimated time remaining:", remainining), "\n")
    invisible(gc())
  }
}

# 5.2. Create bulk data for training:
bulk_options = expand.grid(c(0, 1), c('no', 'limit_1000'), c(110, 1100), c('YES', 'NO'),
                           c(100, 500), c(10000, 20000), c('YES', 'NO'), stringsAsFactors=F)
rownames(bulk_options) = apply(expand.grid(c('A1', 'A2'), c('B1', 'B2'), c('C1', 'C2'), c('D1', 'D2'),
                                           c('E1', 'E2'), c('F1', 'F2'), c('G1', 'G2'), stringsAsFactors=F), 1, paste, collapse='')
colnames(bulk_options) = c('load_min.counts', 'simulation', 'simul_subset.cells', 'simul_proportional',
                           'bulk_n.cells', 'bulk_num.bulk.samples', 'bulk_balanced.type.cells')
# For simulation 'no', arguments 'simul_subset.cells' and 'simul_proportional' are irrelevant:
bulk_options = rbind(bulk_options[bulk_options$simulation!='no',],
                     bulk_options[rownames(unique(bulk_options[bulk_options$simulation=='no',c(1,5,6,7)])),])
# Arguments of bulk_n.cells=500 and bulk_num.bulk.samples=20000 will be removed. Too much computational burden:
bulk_options = bulk_options[!(bulk_options$bulk_n.cells==500 & bulk_options$bulk_num.bulk.samples==20000),]
# ..
n_params_bulk = dim(bulk_options)[1]
pb_bulk = txtProgressBar(min=0, max=n_params_bulk * 2, style=3, width=200, char='=')
init_bulk = numeric(n_params_bulk * 1)
end_bulk = numeric(n_params_bulk * 1)
cr=1
for(ref in c('uncorrected')){#, 'corrected')){
  for(combo_i in 1:n_params_bulk){
    if(cr==1) ti = combo_i
    if(cr==2) ti = n_params + combo_i
    init_bulk[ti] <- Sys.time()
    combo_name = rownames(bulk_options)[combo_i]
    message(paste(combo_name, '|', ti, '/', length(init_bulk) * 1, sep=' '))
    combo = bulk_options[combo_name, ]
    
    ###
    e_combo = substr(combo_name, 9, 10)
    if(e_combo=='E1') next
    ###
    
    bulk_balanced.type.cells = FALSE
    if(combo$bulk_balanced.type.cells=='YES') bulk_balanced.type.cells = TRUE
    
    if(combo$simulation == 'limit_1000'){
      simul_combo = substr(combo_name, 1, 8)
      DDLSorter_obj = readRDS(paste('./Results/1_Parameter_Optimization/train/DigitalDLSorter_utils/simulate_cells/',
                                    ref, '_', simul_combo, '.Rdata', sep=''))
    }
    else{
      DDLSorter_obj = prepare_singlecellexperiment_scReferences(references[[ref]], 'Deconv_cellTypes', 10,
                                                                combo$load_min.counts, 'DDLS')
    }
    invisible(gc())
    # Calculate minimum and maximum proportions of cell-types by sample in reference to create the minmaxprobs argument
    cts = unique(references[[ref]]$metadata[,'Deconv_cellTypes'])
    min_props = rep(1, length(cts))
    max_props = rep(0, length(cts))
    samps = unique(references[[ref]]$metadata[,'Sample'])
    for(samp in samps){
      tab_cts = table(references[[ref]]$metadata[,'Deconv_cellTypes'][references[[ref]]$metadata[,'Sample']==samp])
      samp_cts = rep(0, length(cts))
      names(samp_cts) = cts
      samp_cts[names(tab_cts)] = tab_cts / sum(tab_cts)
      min_props[samp_cts < min_props] = samp_cts[samp_cts < min_props]
      max_props[samp_cts > max_props] = samp_cts[samp_cts > max_props]
    }
    min_props[min_props==0] = 0.01
    minmaxprobs = data.frame(Cell_Type=cts, from=min_props*100, to=max_props*100)
    invisible(gc())
    # Generate bulk:
    DDLSorter_obj = generate_bulk(DDLSorter_obj, minmaxprobs, combo$bulk_num.bulk.samples, combo$bulk_n.cells,
                                                        bulk_balanced.type.cells, 6)
    saveRDS(DDLSorter_obj, paste('./Results/1_Parameter_Optimization/train/DigitalDLSorter_utils/bulk_DDLS/', ref, '_', combo_name, '.Rdata', sep=''))
    remove(DDLSorter_obj)
    invisible(gc())
    
    message('')
    end_bulk[ti] = Sys.time()
    setTxtProgressBar(pb_bulk, ti)
    time = round(lubridate::seconds_to_period(sum(end_bulk - init_bulk)), 0)
    est = n_params_bulk * 1 * (mean(end_bulk[end_bulk != 0] - init_bulk[init_bulk != 0])) - sum(end_bulk - init_bulk)
    remainining <- round(lubridate::seconds_to_period(est), 0)
    cat(paste(" // Execution time:", time, " // Estimated time remaining:", remainining), "\n")
    invisible(gc())
  }
}

# 5.3. Training:
train_rmses = c()
train_rmses_samp = c()
options = c()
models = c()

n_params = dim(methods_parameters$DigitalDLSorter)[1]
pb = txtProgressBar(min=0, max=n_params * 3, style=3, width=200, char='=')
init = numeric(n_params * 3)
end = numeric(n_params * 3)
cr = 1
best_model = list()
for(correction in c('corrected.Before', 'corrected.After', 'uncorrected')){
  for(combo_i in 1:n_params){
    if(cr==1) ti = combo_i
    if(cr==2) ti = n_params + combo_i
    if(cr==3) ti = (n_params*2) + combo_i
    init[ti] <- Sys.time()
    combo_name = rownames(methods_parameters$DigitalDLSorter)[combo_i]
    options = c(options, paste(correction, combo_name, sep='_'))
    message(paste(combo_name, '|', ti, '/', length(init), sep=' '))
    combo = methods_parameters$DigitalDLSorter[combo_name, ]
    
    predict_normalize = FALSE
    if(combo$predict_normalize=='YES') predict_normalize = TRUE
    
    obj_name = substr(combo_name, 1, 14)
    ref = 'uncorrected'
    if(correction=='corrected.Before') ref = corrected
    DDLSorter_obj = readRDS(paste('./Results/1_Parameter_Optimization/train/DigitalDLSorter_utils/bulk_DDLS/',
                                  ref, '_', obj_name, '.Rdata', sep=''))

    # Training:
    DDLSorter_obj = train_DDLSorter_model(DDLSorter_obj, combo$train_combine, 128, 10, 6, FALSE)
    invisible(gc())
    # Predicting:
    res = deconv_DigitalDLSorter(DDLSorter_obj, train, 128, predict_normalize, 'deconv')
    if(correction=='corrected.After') res = correct_fractions_mRNABias(res, rna_correction_vec)
    invisible(gc())
    
    rmses_samp = c()
    for(samp in colnames(res)) rmses_samp = c(rmses_samp, Metrics::rmse(unlist(train$metadata[cell_types, samp]),
                                                                        unlist(res[cell_types, samp])))
    train_rmses_samp = cbind(train_rmses_samp, rmses_samp)
    train_rmses = c(train_rmses, mean(rmses_samp))
    
    message('Storing best model')
    if(ti == 1){
      best_model$rmse = mean(rmses_samp)
      best_model$model = DDLSorter_obj
    }
    else{
      if(best_model$rmse > mean(rmses_samp)){
        best_model$rmse = mean(rmses_samp)
        best_model$model = DDLSorter_obj
      }
    }
    remove(DDLSorter_obj)
    invisible(gc())
    
    message('')
    end[ti] = Sys.time()
    setTxtProgressBar(pb, ti)
    time = round(lubridate::seconds_to_period(sum(end - init)), 0)
    est = n_params * 3 * (mean(end[end != 0] - init[init != 0])) - sum(end - init)
    remainining <- round(lubridate::seconds_to_period(est), 0)
    cat(paste(" // Execution time:", time, " // Estimated time remaining:", remainining), "\n")
  }
  remove(reference)
  invisible(gc())
  cr = cr+1
}
names(train_rmses) = options
colnames(train_rmses_samp) = options
rownames(train_rmses_samp) = colnames(res)
names(models) = models

# 5.4. Save train results in csv:
full_train_results = rbind(train_rmses_samp, train_rmses, time)
nrows = dim(full_train_results)[1]
rownames(full_train_results)[(nrows-1):nrows] = c('MEAN', 'TIME')
write.csv(full_train_results, paste('./Results/1_Parameter_Optimization/train/DigitalDLSorter.csv', sep=''))



# 6. DWLS
train_rmses = c()
train_rmses_samp = c()
options = c()

# 6.1. Training:
n_params = dim(methods_parameters$DWLS)[1]
pb = txtProgressBar(min=0, max=n_params * 3, style=3, width=n_params * 3, char='=')
init = numeric(n_params * 2)
end = numeric(n_params * 2)
cr = 1
for(correction in c('corrected.After', 'uncorrected')){
  for(combo_i in 1:n_params){
    if(cr==1) ti = combo_i
    if(cr==2) ti = n_params + combo_i
    init[ti] <- Sys.time()
    combo_name = rownames(methods_parameters$DWLS)[combo_i]
    options = c(options, paste(correction, combo_name, sep='_'))
    message(paste(combo_name, '|', ti, '/', length(init), sep=' '))
    combo = methods_parameters$DWLS[combo_name, ]
    res = Estimate.Proportions(train, references[['uncorrected']], 'Deconv_cellTypes', 'DWLS', combo$diff.cutoff,
                               combo$pval.cutoff, rna.bias=FALSE, bias.vec=NULL)
    if(correction=='corrected.After') res = correct_fractions_mRNABias(res, rna_correction_vec)
    rmses_samp = c()
    for(samp in colnames(res)) rmses_samp = c(rmses_samp, Metrics::rmse(unlist(train$metadata[cell_types, samp]),
                                                                        unlist(res[cell_types, samp])))
    train_rmses_samp = cbind(train_rmses_samp, rmses_samp)
    train_rmses = c(train_rmses, mean(rmses_samp))
    message('')
    end[ti] = Sys.time()
    setTxtProgressBar(pb, ti)
    time = round(lubridate::seconds_to_period(sum(end - init)), 0)
    est = n_params * 3 * (mean(end[end != 0] - init[init != 0])) - sum(end - init)
    remainining <- round(lubridate::seconds_to_period(est), 0)
    cat(paste(" // Execution time:", time, " // Estimated time remaining:", remainining), "\n")
    invisible(gc())
  }
  invisible(gc())
  cr = cr+1
}
names(train_rmses) = options
time = end - init
names(time) = options
colnames(train_rmses_samp) = options
rownames(train_rmses_samp) = colnames(res)

# 6.2. Save train results in csv:
full_train_results = rbind(train_rmses_samp, train_rmses, time)
nrows = dim(full_train_results)[1]
rownames(full_train_results)[(nrows-1):nrows] = c('MEAN', 'TIME')
write.csv(full_train_results, './Results/1_Parameter_Optimization/train/DWLS.csv')



# 7. MOMF
train_rmses = c()
train_rmses_samp = c()
options = c()

# 7.1. Training:
n_params = dim(methods_parameters$MOMF)[1]
pb = txtProgressBar(min=0, max=n_params * 3, style=3, width=n_params * 3, char='=')
init = numeric(n_params * 3)
end = numeric(n_params * 3)
cr = 1
for(correction in c('corrected.Before', 'corrected.After', 'uncorrected')){
  if(correction=='corrected.Before') reference = references[['corrected']]
  else reference = references[['uncorrected']]
  for(combo_i in 1:n_params){
    if(cr==1) ti = combo_i
    if(cr==2) ti = n_params + combo_i
    if(cr==3) ti = (n_params*2) + combo_i
    init[ti] <- Sys.time()
    combo_name = rownames(methods_parameters$MOMF)[combo_i]
    options = c(options, paste(correction, combo_name, sep='_'))
    message(paste(combo_name, '|', ti, '/', length(init), sep=' '))
    combo = methods_parameters$MOMF[combo_name, ]
    res = Estimate.Proportions(train, reference, 'Deconv_cellTypes', 'MOMF', combo$method, combo$rho, combo$num.iter, FALSE, NULL)
    if(correction=='corrected.After') res = correct_fractions_mRNABias(res, rna_correction_vec)
    rmses_samp = c()
    for(samp in colnames(res)) rmses_samp = c(rmses_samp, Metrics::rmse(unlist(train$metadata[cell_types, samp]),
                                                                        unlist(res[cell_types, samp])))
    train_rmses_samp = cbind(train_rmses_samp, rmses_samp)
    train_rmses = c(train_rmses, mean(rmses_samp))
    message('')
    end[ti] = Sys.time()
    setTxtProgressBar(pb, ti)
    time = round(lubridate::seconds_to_period(sum(end - init)), 0)
    est = n_params * 3 * (mean(end[end != 0] - init[init != 0])) - sum(end - init)
    remainining <- round(lubridate::seconds_to_period(est), 0)
    cat(paste(" // Execution time:", time, " // Estimated time remaining:", remainining), "\n")
    invisible(gc())
  }
  remove(reference)
  invisible(gc())
  cr = cr+1
}
names(train_rmses) = options
time = end - init
names(time) = options
colnames(train_rmses_samp) = options
rownames(train_rmses_samp) = colnames(res)

# 7.2. Save train results in csv:
full_train_results = rbind(train_rmses_samp, train_rmses, time)
nrows = dim(full_train_results)[1]
rownames(full_train_results)[(nrows-1):nrows] = c('MEAN', 'TIME')
write.csv(full_train_results, './Results/1_Parameter_Optimization/train/MOMF.csv')



# 8. MuSiC_woGrouping
train_rmses = c()
train_rmses_samp = c()
options = c()

# 8.1. Training:
n_params = dim(methods_parameters$MuSiC_woGrouping)[1]
pb = txtProgressBar(min=0, max=n_params * 2, style=3, width=n_params * 2, char='=')
init = numeric(n_params)
end = numeric(n_params)
reference = references[['uncorrected']]
for(combo_i in 1:n_params){
  init[combo_i] <- Sys.time()
  combo_name = rownames(methods_parameters$MuSiC_woGrouping)[combo_i]
  options = c(options, paste(combo_name, sep='_'))
  message(paste(combo_name, '|', combo_i, '/', length(init), sep=' '))
  combo = methods_parameters$MuSiC_woGrouping[combo_name, ]
  gene_markers = NULL
  if(combo$wo.gene.markers=='YES') gene_markers = markers
  method = 'Est.prop.weighted'
  if(combo$wo.method=='all genes') method = 'Est.prop.allgene'
  res = Estimate.Proportions(train, reference, 'Deconv_cellTypes', 'MuSiC', music.method='without.Grouping', n.iterations=combo$n.iterations,
                             center=combo$center, normalize=combo$normalise, wo.method=method, wo.gene.markers=gene_markers,
                             wo.bias.vec=rna_correction_vec)
  rmses_samp = c()
  for(samp in colnames(res)) rmses_samp = c(rmses_samp, Metrics::rmse(unlist(train$metadata[cell_types, samp]),
                                                                      unlist(res[cell_types, samp])))
  train_rmses_samp = cbind(train_rmses_samp, rmses_samp)
  train_rmses = c(train_rmses, mean(rmses_samp))
  message('')
  end[combo_i] = Sys.time()
  setTxtProgressBar(pb, combo_i)
  time = round(lubridate::seconds_to_period(sum(end - init)), 0)
  est = n_params * (mean(end[end != 0] - init[init != 0])) - sum(end - init)
  remainining <- round(lubridate::seconds_to_period(est), 0)
  cat(paste(" // Execution time:", time, " // Estimated time remaining:", remainining), "\n")
  invisible(gc())
}
names(train_rmses) = options
time = end - init
names(time) = options
colnames(train_rmses_samp) = options
rownames(train_rmses_samp) = colnames(res)

# 8.2. Save train results in csv:
full_train_results = rbind(train_rmses_samp, train_rmses, time)
nrows = dim(full_train_results)[1]
rownames(full_train_results)[(nrows-1):nrows] = c('MEAN', 'TIME')
write.csv(full_train_results, './Results/1_Parameter_Optimization/train/MuSiC_woGrouping.csv')



# 9. MuSiC_wGrouping
train_rmses = c()
train_rmses_samp = c()
options = c()

# 9.1. Training:
n_params = dim(methods_parameters$MuSiC_wGrouping)[1]
pb = txtProgressBar(min=0, max=n_params * 3, style=3, width=n_params * 3, char='=')
init = numeric(n_params * 3)
end = numeric(n_params * 3)
cr = 1
clusters = list(C1=c("Cancer cells", "Stromal cells"), C2=c("Anti-Inflammatory macro/mono", "Pro-Inflammatory macro/mono", "Other cells"),
                C3=c("Bcells", "Regulatory CD4 Tcells", "CD4 Tcells", "Proliferative Tcells", "NK cells", "CD8 Tcells"))
# Run algorithm:
for(correction in c('corrected.Before', 'corrected.After', 'uncorrected')){
  if(correction=='corrected.Before') reference = references[['corrected']]
  else reference = references[['uncorrected']]
  scReference = prepare_expressionSet_scReferences(reference, 'Deconv_cellTypes')
  bulk_data_expressionSet = Biobase::ExpressionSet(data.matrix(train$data))
  exprObj = Seurat::CreateSeuratObject(as.data.frame(exprs(scReference)), project = "DE")
  exprObj = Seurat::SetIdent(exprObj, value=as.vector(scReference$cellType))
  exprObj = Seurat::SCTransform(exprObj)
  invisible(gc())
  # Calculate markers for clusters:
  group_markers=list()
  for(cl in names(clusters)){
    message(paste('- Cluster:', cl, sep=' '))
    if(length(clusters[[cl]])>1){
      markers = c()
      for(cell_type in clusters[[cl]]){
        message('--', cell_type, ' ')
        other_cells = clusters[[cl]][-match(cell_type, clusters[[cl]])]
        cell_type_markers = Seurat::FindMarkers(exprObj, ident.1=cell_type, ident.2=other_cells, slot='counts',
                                                test.use='wilcox', logfc.threshold=0.8, min.pct=0.3, only.pos=TRUE)
        markers = c(markers, rownames(cell_type_markers)[cell_type_markers$p_val_adj<0.01])
      }
      group_markers[[cl]] = unique(markers)
    }
    else{
      group_markers[cl] = list(NULL)
      message('Only one cell-type in this cluster')
    }
  }
  for(combo_i in 1:n_params){
    if(cr==1) ti = combo_i
    if(cr==2) ti = n_params + combo_i
    if(cr==3) ti = (n_params*2) + combo_i
    init[ti] <- Sys.time()
    combo_name = rownames(methods_parameters$MuSiC_wGrouping)[combo_i]
    options = c(options, paste(correction, combo_name, sep='_'))
    message(paste(combo_name, '|', ti, '/', length(init), sep=' '))
    combo = methods_parameters$MuSiC_wGrouping[combo_name, ]
    # Estimate proportions:
    clusterType = as.character(scReference$cellType)
    for(cl in 1:length(clusters)){
      clusterType[clusterType%in%clusters[[cl]]] = names(clusters)[cl]
    }
    pData(scReference)$clusterType = factor(clusterType, levels=names(clusters))
    music_proportions = MuSiC::music_prop.cluster(bulk.eset=bulk_data_expressionSet, sc.eset=scReference,
                                                  group.markers=group_markers, groups='clusterType',
                                                  clusters='cellType', samples='sampleID',
                                                  clusters.type=clusters,
                                                  verbose=T, iter.max=combo$n.iterations, centered=combo$center, normalize=combo$normalise)
    res = t(music_proportions$Est.prop.weighted.cluster)
    if(correction=='corrected.After') res = correct_fractions_mRNABias(res, rna_correction_vec)
    rmses_samp = c()
    for(samp in colnames(res)) rmses_samp = c(rmses_samp, Metrics::rmse(unlist(train$metadata[cell_types, samp]),
                                                                        unlist(res[cell_types, samp])))
    train_rmses_samp = cbind(train_rmses_samp, rmses_samp)
    train_rmses = c(train_rmses, mean(rmses_samp))
    message('')
    end[ti] = Sys.time()
    setTxtProgressBar(pb, ti)
    time = round(lubridate::seconds_to_period(sum(end - init)), 0)
    est = n_params * 3 * (mean(end[end != 0] - init[init != 0])) - sum(end - init)
    remainining <- round(lubridate::seconds_to_period(est), 0)
    cat(paste(" // Execution time:", time, " // Estimated time remaining:", remainining), "\n")
    invisible(gc())
  }
  remove(reference, exprObj)
  invisible(gc())
  cr = cr+1
}
names(train_rmses) = options
time = end - init
names(time) = options
colnames(train_rmses_samp) = options
rownames(train_rmses_samp) = colnames(res)

# 9.2. Save train results in csv:
full_train_results = rbind(train_rmses_samp, train_rmses, time)
nrows = dim(full_train_results)[1]
rownames(full_train_results)[(nrows-1):nrows] = c('MEAN', 'TIME')
write.csv(full_train_results, './Results/1_Parameter_Optimization/train/MuSiC_wGrouping.csv')



# 10. Scaden
train_rmses = c()
train_rmses_samp = c()
options = c()

# 10.1. Training:
n_params = dim(methods_parameters$Scaden)[1]
pb = txtProgressBar(min=0, max=n_params * 3, style=3, width=n_params * 3, char='=')
init = numeric(n_params * 3)
end = numeric(n_params * 3)
cr = 1
cell_types_scaden = cell_types
names(cell_types_scaden) = cell_types
cell_types_scaden[3:4] = c('Anti Inflammatory macro mono', 'Pro Inflammatory macro mono')
rna_correction_vec_scaden = rna_correction_vec
names(rna_correction_vec_scaden) = cell_types_scaden[match(names(rna_correction_vec), names(cell_types_scaden))]
for(correction in c('corrected.Before', 'corrected.After', 'uncorrected')){
  if(correction=='corrected.Before') train_file = scaden_train_files[['corrected']]
  else train_file = scaden_train_files[['uncorrected']]
  for(combo_i in 1:n_params){
    if(cr==1) ti = combo_i
    if(cr==2) ti = n_params + combo_i
    if(cr==3) ti = (n_params*2) + combo_i
    init[ti] <- Sys.time()
    combo_name = rownames(methods_parameters$Scaden)[combo_i]
    options = c(options, paste(correction, combo_name, sep='_'))
    message(paste(combo_name, '|', ti, '/', length(init), sep=' '))
    combo = methods_parameters$Scaden[combo_name, ]
    res = Estimate.Proportions(train_bulk_file_scaden, train_file, cellType.var, 'Scaden',
                               dir.results=tempdir(), datasets='', min.expression=combo$min.expression, batch.size=as.integer(combo$batch.size),
                               learning.rate=combo$learning.rate, n.steps=as.integer(combo$n.steps))
    if(correction=='corrected.After') res = correct_fractions_mRNABias(res, rna_correction_vec_scaden)
    rmses_samp = c()
    for(samp in colnames(res)) rmses_samp = c(rmses_samp, Metrics::rmse(unlist(train$metadata[cell_types, samp]),
                                                                        unlist(res[cell_types_scaden, samp])))
    train_rmses_samp = cbind(train_rmses_samp, rmses_samp)
    train_rmses = c(train_rmses, mean(rmses_samp))
    message('')
    end[ti] = Sys.time()
    setTxtProgressBar(pb, ti)
    time = round(lubridate::seconds_to_period(sum(end - init)), 0)
    est = n_params * 3 * (mean(end[end != 0] - init[init != 0])) - sum(end - init)
    remainining <- round(lubridate::seconds_to_period(est), 0)
    cat(paste(" // Execution time:", time, " // Estimated time remaining:", remainining), "\n")
    invisible(gc())
  }
  invisible(gc())
  cr = cr+1
}
names(train_rmses) = options
time = end - init
names(time) = options
colnames(train_rmses_samp) = options
rownames(train_rmses_samp) = colnames(res)

# 10.2. Save train results in csv:
full_train_results = rbind(train_rmses_samp, train_rmses, time)
nrows = dim(full_train_results)[1]
rownames(full_train_results)[(nrows-1):nrows] = c('MEAN', 'TIME')
write.csv(full_train_results, './Results/1_Parameter_Optimization/train/Scaden.csv')



# 11. SCDC
train_rmses = c()
train_rmses_samp = c()
options = c()

# 11.1. Training:
n_params = dim(methods_parameters$SCDC)[1]
pb = txtProgressBar(min=0, max=n_params * 3, style=3, width=n_params * 3, char='=')
init = numeric(n_params * 3)
end = numeric(n_params * 3)
cr = 1
for(correction in c('corrected.Before', 'corrected.After', 'uncorrected')){
  if(correction=='corrected.Before') reference = references[['corrected']]
  else reference = references[['uncorrected']]
  for(combo_i in 1:n_params){
    if(cr==1) ti = combo_i
    if(cr==2) ti = n_params + combo_i
    if(cr==3) ti = (n_params*2) + combo_i
    init[ti] <- Sys.time()
    combo_name = rownames(methods_parameters$SCDC)[combo_i]
    options = c(options, paste(correction, combo_name, sep='_'))
    message(paste(combo_name, '|', ti, '/', length(init), sep=' '))
    combo = methods_parameters$SCDC[combo_name, ]
    res = Estimate.Proportions(train, reference, 'Deconv_cellTypes', 'SCDC', scdc.method='multiple', n.iterations=combo$m.n.iterations,
                               m.use.gridSearch=FALSE, m.search.length=combo$m.search.length, m.dataset.var='Dataset', m.weights.method=NULL)
    for(method in names(res)){
      if(correction=='corrected.After') res[[method]] = correct_fractions_mRNABias(res[[method]], rna_correction_vec)
      rmses_samp = c()
      cts = cell_types[!is.na(match(cell_types, rownames(res[[method]])))]
      for(samp in colnames(res[[method]]))
        rmses_samp = c(rmses_samp, Metrics::rmse(unlist(train$metadata[cts, samp]), unlist(res[[method]][cts, samp])))
      train_rmses_samp[[method]] = cbind(train_rmses_samp[[method]], rmses_samp)
      train_rmses[[method]] = c(train_rmses[[method]], mean(rmses_samp))
    }
    message('')
    end[ti] = Sys.time()
    setTxtProgressBar(pb, ti)
    time = round(lubridate::seconds_to_period(sum(end - init)), 0)
    est = n_params * 3 * (mean(end[end != 0] - init[init != 0])) - sum(end - init)
    remainining <- round(lubridate::seconds_to_period(est), 0)
    cat(paste(" // Execution time:", time, " // Estimated time remaining:", remainining), "\n")
    invisible(gc())
  }
  remove(reference)
  invisible(gc())
  cr = cr+1
}
for(method in names(train_rmses)) names(train_rmses[[method]]) = options
for(method in names(train_rmses_samp)) colnames(train_rmses_samp[[method]]) = options
for(method in names(train_rmses_samp)) rownames(train_rmses_samp[[method]]) = colnames(res[[1]])

# 11.2. Save train results in csv:
full_train_results = list()
for(method in names(train_rmses_samp)){
  full_train_results[[method]] = rbind(train_rmses_samp[[method]], train_rmses[[method]])
  nrows = dim(full_train_results[[method]])[1]
  rownames(full_train_results[[method]])[nrows] = 'MEAN'
}
write.csv(full_train_results, paste('./Results/1_Parameter_Optimization/train/SCDC.csv', sep=''))
