code_dir = './code/R'
for(file in list.files(code_dir, full.names=TRUE)) source(file)

methods = c('AutoGeneS_linear', 'AutoGeneS_nnls', 'AutoGeneS_nusvr', 'BisqueRNA', 'BSeqSC', 'CIBERSORTx', 'DigitalDLSorter', 'DWLS', 'MOMF',
            'MuSiC_woGrouping', 'MuSiC_wGrouping', 'Scaden', 'SCDC')
cell_types = c('Cancer cells', 'Stromal cells', 'Macro/mono Lineage', 'Bcells', 'CD4 Tcells',
               'Regulatory CD4 Tcells', 'CD8 Tcells','Proliferative Tcells', 'NK cells')

methods_colours = c('#808000', '#5c4033', '#f2d2bd', '#c4961a', '#b30000', '#d16103', '#ffdb6d',
                    '#c3d7a4','#52854c', '#4e84c4','#293352','#999999', '#cc79a7')
names(methods_colours) = methods





# --------------------
# --- READ RESULTS ---
# --------------------

results_dir = './Results/proportions'

# Read result files:
results = list()
for(type in c('uncorrected', 'corrected_before', 'corrected_after', '')){
  x = type
  if(type==''){
    results$other_models = list()
    x = 'other_models'
  }
  else results[[type]] = list()
  results[[x]]$woProl = list()
  results[[x]]$wProl = list()
  
  res_files = grep('[.]csv', list.files(paste(results_dir, type, sep='/')), value=TRUE)
  names(res_files) = rep('wProl', length(res_files))
  names(res_files)[grep('woProl_', res_files)] = 'woProl'
  
  for(file_idx in 1:length(res_files)){
    if(type=='') res = read.csv(paste(results_dir, res_files[file_idx], sep='/'), row.names=1)
    else res = read.csv(paste(results_dir, type, res_files[file_idx], sep='/'), row.names=1)
    model_name =  gsub('[.]csv', '', gsub('woProl_', '', res_files[file_idx]))
    results[[x]][[names(res_files)[file_idx]]][[model_name]] = res
  }
}

# Exclude nova samples:
samples_to_keep = c("NIC12", "NIC13", "NIC15", "NIC16", "NIC17", "NIC18", "NIC19", "NIC20", "NIC21", "NIC22", "NIC23",
                    "NIC24", "NIC25", "NIC27", "NIC29", "NIC3", "NIC4", "NIC5", "NIC6", "NIC7", "S52")
for(correction in names(results)){
  for(prol_status in names(results[[correction]])){
    for(method in names(results[[correction]][[prol_status]])){
      results[[correction]][[prol_status]][[method]] = results[[correction]][[prol_status]][[method]][,samples_to_keep]
    }
  }
}





# -------------------------------------
# --- READ GROUND TRUTH PROPORTIONS ---
# -------------------------------------

ground_truth_file_wOther = './Data/bulk/CRC/NICs/data/cell_proportions.csv'
ground_truth_file_woOther = './Data/bulk/CRC/NICs/data/cell_proportions_woOthers.csv'

ground_truth = read.csv(ground_truth_file_wOther, row.names=1)[, samples_to_keep]
ground_truth_woOther = read.csv(ground_truth_file_woOther, row.names=1)[, samples_to_keep]





# ------------------------------------------------------------
# --- CALCULATE RMSE AND CORRELATION OF MODELS PREDICTIONS ---
# ------------------------------------------------------------


# 0. Create matrices where average RMSEs and correlations will be stored and metadata:
rmses = matrix(data=NA, nrow=length(methods), ncol=12, dimnames=list(methods, c('uncorrected_gd', 'correctedBefore_gd', 'correctedAfter_gd',
                                                                                'uncorrected_woProl_gd', 'correctedBefore_woProl_gd',
                                                                                'correctedAfter_woProl_gd', 'uncorrected_gdwoOthers',
                                                                                'correctedBefore_gdwoOthers', 'correctedAfter_gdwoOthers',
                                                                                'uncorrected_woProl_gdwoOthers', 'correctedBefore_woProl_gdwoOthers',
                                                                                'correctedAfter_woProl_gdwoOthers')))
corrs = matrix(data=NA, nrow=length(methods), ncol=12, dimnames=list(methods, c('uncorrected_gd', 'correctedBefore_gd', 'correctedAfter_gd',
                                                                                'uncorrected_woProl_gd', 'correctedBefore_woProl_gd',
                                                                                'correctedAfter_woProl_gd', 'uncorrected_gdwoOthers',
                                                                                'correctedBefore_gdwoOthers', 'correctedAfter_gdwoOthers',
                                                                                'uncorrected_woProl_gdwoOthers', 'correctedBefore_woProl_gdwoOthers',
                                                                                'correctedAfter_woProl_gdwoOthers')))
metadata = data.frame(reference = c('uncorrected', 'correctedBefore', 'correctedAfter', 'uncorrected', 'correctedBefore', 'correctedAfter'),
                      proliferative = c(rep('wProl', 3), rep('woProl', 3), rep('wProl', 3), rep('woProl', 3)),
                      ground_truth = c(rep('gdwOthers', 6), rep('woOthers', 6)),
                      row.names=c('uncorrected_gd', 'correctedBefore_gd', 'correctedAfter_gd',
                                  'uncorrected_woProl_gd', 'correctedBefore_woProl_gd',
                                  'correctedAfter_woProl_gd', 'uncorrected_gdwoOthers',
                                  'correctedBefore_gdwoOthers', 'correctedAfter_gdwoOthers',
                                  'uncorrected_woProl_gdwoOthers', 'correctedBefore_woProl_gdwoOthers',
                                  'correctedAfter_woProl_gdwoOthers'))

# 1. Calculate RMSEs:
samps_rmses = c()
samps_names = c()
samps_models = c()
samps_cols = c()
cts_rmses = c()
cts_names = c()
cts_models = c()
cts_cols = c()

for(correction in names(results)){
  x = 'uncorrected'
  if(correction=='corrected_before') x = 'correctedBefore'
  else if(correction=='corrected_after') x = 'correctedAfter'
  
  for(prol_type in names(results[[correction]])){
    if(prol_type=='woProl') y = paste(x, prol_type, sep='_')
    else y = x
    cts = cell_types
    if(prol_type == 'woProl') cts = cell_types[-8]
    
    for(model_name in names(results[[correction]][[prol_type]])){
      model_results = results[[correction]][[prol_type]][[model_name]]
      
      rmse_model = rmse_dataset(model_results, ground_truth, cts)
      rmses[model_name, paste(y, 'gd', sep='_')] = rmse_model$dataset_average
      samps_rmses = c(samps_rmses, rmse_model$samples)
      samps_names = c(samps_names, names(rmse_model$samples))
      samps_models = c(samps_models, rep(model_name, length(rmse_model$samples)))
      samps_cols = c(samps_cols, rep(paste(y, 'gd', sep='_'), length(rmse_model$samples)))
      cts_rmses = c(cts_rmses, rmse_model$cell_types)
      cts_names = c(cts_names, names(rmse_model$cell_types))
      cts_models = c(cts_models, rep(model_name, length(rmse_model$cell_types)))
      cts_cols = c(cts_cols, rep(paste(y, 'gd', sep='_'), length(rmse_model$cell_types)))
      
      rmse_model2 = rmse_dataset(model_results, ground_truth_woOther, cts)
      rmses[model_name, paste(y, 'gdwoOthers', sep='_')] = rmse_model2$dataset_average
      samps_rmses = c(samps_rmses, rmse_model2$samples)
      samps_names = c(samps_names, names(rmse_model2$samples))
      samps_models = c(samps_models, rep(model_name, length(rmse_model2$samples)))
      samps_cols = c(samps_cols, rep(paste(y, 'gdwoOthers', sep='_'), length(rmse_model2$samples)))
      cts_rmses = c(cts_rmses, rmse_model2$cell_types)
      cts_names = c(cts_names, names(rmse_model2$cell_types))
      cts_models = c(cts_models, rep(model_name, length(rmse_model2$cell_types)))
      cts_cols = c(cts_cols, rep(paste(y, 'gdwoOthers', sep='_'), length(rmse_model2$cell_types)))
    }
  }
}
rmses_samples_df = data.frame(values=samps_rmses, samples=samps_names, models=samps_models, cols=samps_cols)
rmses_cellTypes_df = data.frame(values=cts_rmses, cell_types=cts_names, models=cts_models, cols=cts_cols)

write.csv(rmses, './Results/metrics/average_rmses.csv')
write.csv(rmses_samples_df, './Results/metrics/samples_rmses.csv')
write.csv(rmses_cellTypes_df, './Results/metrics/cellTypes_rmses.csv')


# 1. Calculate correlations:
samps_corrs = c()
samps_names = c()
samps_models = c()
samps_cols = c()
cts_corrs = c()
cts_names = c()
cts_models = c()
cts_cols = c()

for(correction in names(results)){
  x = 'uncorrected'
  if(correction=='corrected_before') x = 'correctedBefore'
  else if(correction=='corrected_after') x = 'correctedAfter'
  
  for(prol_type in names(results[[correction]])){
    if(prol_type=='woProl') y = paste(x, prol_type, sep='_')
    else y = x
    cts = cell_types
    if(prol_type == 'woProl') cts = cell_types[-8]
    
    for(model_name in names(results[[correction]][[prol_type]])){
      model_results = results[[correction]][[prol_type]][[model_name]]
      
      corr_model = correlation_dataset(model_results, ground_truth, cts)
      corrs[model_name, paste(y, 'gd', sep='_')] = corr_model$dataset_average
      samps_corrs = c(samps_corrs, corr_model$samples)
      samps_names = c(samps_names, names(corr_model$samples))
      samps_models = c(samps_models, rep(model_name, length(corr_model$samples)))
      samps_cols = c(samps_cols, rep(paste(y, 'gd', sep='_'), length(corr_model$samples)))
      cts_corrs = c(cts_corrs, corr_model$cell_types)
      cts_names = c(cts_names, names(corr_model$cell_types))
      cts_models = c(cts_models, rep(model_name, length(corr_model$cell_types)))
      cts_cols = c(cts_cols, rep(paste(y, 'gd', sep='_'), length(corr_model$cell_types)))
      
      corr_model2 = correlation_dataset(model_results, ground_truth_woOther, cts)
      corrs[model_name, paste(y, 'gdwoOthers', sep='_')] = corr_model2$dataset_average
      samps_corrs = c(samps_corrs, corr_model2$samples)
      samps_names = c(samps_names, names(corr_model2$samples))
      samps_models = c(samps_models, rep(model_name, length(corr_model2$samples)))
      samps_cols = c(samps_cols, rep(paste(y, 'gdwoOthers', sep='_'), length(corr_model2$samples)))
      cts_corrs = c(cts_corrs, corr_model2$cell_types)
      cts_names = c(cts_names, names(corr_model2$cell_types))
      cts_models = c(cts_models, rep(model_name, length(corr_model2$cell_types)))
      cts_cols = c(cts_cols, rep(paste(y, 'gdwoOthers', sep='_'), length(corr_model2$cell_types)))
    }
  }
}
correlations_samples_df = data.frame(values=samps_corrs, samples=samps_names, models=samps_models, cols=samps_cols)
correlations_cellTypes_df = data.frame(values=cts_corrs, cell_types=cts_names, models=cts_models, cols=cts_cols)

write.csv(corrs, './Results/metrics/average_correlations.csv')
write.csv(correlations_samples_df, './Results/metrics/samples_correlations.csv')
write.csv(correlations_cellTypes_df, './Results/metrics/cellTypes_correlations.csv')

