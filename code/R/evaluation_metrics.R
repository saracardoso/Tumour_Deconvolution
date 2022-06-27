
# ------------
# --- RMSE ---
# ------------

rmse_sample = function(predicted_vec, ground_truth_vec, cell_types){
  n_cts = length(cell_types)
  if((sum(cell_types %in% names(predicted_vec)) + sum(cell_types %in% names(ground_truth_vec))) != (n_cts*2)) stop('Incorrect cell_types given')
  
  res = Metrics::rmse(ground_truth_vec[cell_types], predicted_vec[cell_types])
  return(res)
}

rmse_dataset = function(predicted_props_df, ground_truth_df, cell_types){
  samps = colnames(predicted_props_df)
  if(sum(!samps%in%colnames(ground_truth_df)) > 0) stop('Predicted and ground-truth samples do not match!')
  n_cts = length(cell_types)
  if((sum(cell_types %in% rownames(predicted_props_df)) + sum(cell_types %in% rownames(ground_truth_df))) != (n_cts*2))
    stop('Incorrect cell_types given')
  
  rmses = list()
  rmses$samples = c()
  for(samp in samps){
    rmses$samples = c(rmses$samples, Metrics::rmse(predicted_props_df[cell_types,samp], ground_truth_df[cell_types, samp]))
  }
  names(rmses$samples) = samps
  
  rmses$cell_types = c()
  for(ct in cell_types){
    rmses$cell_types = c(rmses$cell_types, Metrics::rmse(as.numeric(predicted_props_df[ct, ]),
                                                         as.numeric(ground_truth_df[ct, colnames(predicted_props_df)])))
  }
  names(rmses$cell_types) = cell_types
  
  rmses$dataset_average = mean(rmses$samples)
  
  return(rmses)
}




# -------------------
# --- Correlation ---
# -------------------

correlation_sample = function(predicted_vec, ground_truth_vec, cell_types){
  n_cts = length(cell_types)
  if((sum(cell_types %in% names(predicted_vec)) + sum(cell_types %in% names(ground_truth_vec))) != (n_cts*2)) stop('Incorrect cell_types given')
  
  res = cor(ground_truth_vec[cell_types], predicted_vec[cell_types], method = 'pearson')
  return(res)
}

correlation_dataset = function(predicted_props_df, ground_truth_df, cell_types){
  samps = colnames(predicted_props_df)
  if(sum(!samps%in%colnames(ground_truth_df)) > 0) stop('Predicted and ground-truth samples do not match!')
  n_cts = length(cell_types)
  if((sum(cell_types %in% rownames(predicted_props_df)) + sum(cell_types %in% rownames(ground_truth_df))) != (n_cts*2))
    stop('Incorrect cell_types given')
  
  correlations = list()
  correlations$samples = c()
  for(samp in samps){
    correlations$samples = c(correlations$samples, cor(predicted_props_df[cell_types,samp], ground_truth_df[cell_types, samp],
                                                       method='pearson'))
  }
  names(correlations$samples) = samps
  
  correlations$cell_types = c()
  for(ct in cell_types){
    correlations$cell_types = c(correlations$cell_types, cor(as.numeric(predicted_props_df[ct, ]),
                                                             as.numeric(ground_truth_df[ct, colnames(predicted_props_df)])))
  }
  names(correlations$cell_types) = cell_types
  
  correlations$dataset_average = mean(correlations$samples)
  
  return(correlations)
}

