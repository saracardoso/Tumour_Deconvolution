
prepare_data_ml = function(seurat_object, n_features=700, assay='SCT', slot='data', cell_type_variable='level1'){
   # Get most variable features:
  seurat_object = Seurat::FindVariableFeatures(seurat_object, nfeatures=n_features, assay=assay)
  
  # Create data to use in training, only with variable features:
  seurat_data = as.matrix(Seurat::GetAssayData(seurat_object, assay=assay, slot=slot)[Seurat::VariableFeatures(seurat_object), ])
  seurat_data = as.data.frame(t(seurat_data))
  seurat_data$class = as.factor(seurat_object@meta.data[rownames(seurat_data), cell_type_variable])
  return(seurat_data)
}


train_models = function(prepared_data, n_folds=3, n_param_opt=5){
  models = list(LDA='lda', PLS='pls', RF='ranger')
  final_results = list()
  n_features = dim(prepared_data)[2]-1
  tictoc::tic.clearlog()
  for(mdl in names(models)){
    message(paste('Training model', mdl))
    tictoc::tic(mdl)
    # Train model:
    train_control = caret::trainControl(method='cv', number=n_folds, summaryFunction=caret::defaultSummary, allowParallel=FALSE)
    if(mdl=='RF') train_result = caret::train(prepared_data[,1:n_features], prepared_data$class, method=models[[mdl]], tuneLength=n_param_opt,
                                                  metric='Accuracy', trControl=train_control, importance='impurity_corrected')
    else train_result = caret::train(prepared_data[,1:n_features], prepared_data$class, method=models[[mdl]], tuneLength=n_param_opt,
                                     metric='Accuracy', trControl=train_control)
    # Organized results:
    final_result = list()
    final_result$var_importance = caret::varImp(train_result)$importance
    final_result$confusion_matrix = caret::confusionMatrix(train_result)
    final_result$optimized_parameters = train_result$bestTune
    final_result$final_model = train_result$finalModel
    final_results[[mdl]] = final_result
    tictoc::toc(log=TRUE, quiet=TRUE)
  }
  log.txt = tictoc::tic.log(format=TRUE)
  message('Time spent running each model:')
  message(paste(unlist(log.txt), collapse='\n'))
  
  return(final_results)
}


predict_full_data = function(train_results, prepared_data){
  
  # Initiate results variables:
  res = list()
  accuracies = c()
  confusion_matrices = list()
  predictions = data.frame(annotated=prepared_data$class, row.names=row.names(prepared_data))
  n_misclassifications = c()
  cell_types = unique(prepared_data$class)
  
  # Run:
  for(mdl in names(train_results)){
    message(paste('Predicting data using model', mdl))
    # Make prediction:
    predict_result = predict(train_results[[mdl]]$final_model, prepared_data[,-dim(prepared_data)[2]])
    if(is.list(predict_result)){
      if(mdl == 'RF') predict_result = predict_result$predictions
      else predict_result = predict_result$class
    }
    
    # Calculate confusion matrix:
    conf_mtx = caret::confusionMatrix(as.factor(predict_result), prepared_data[,'class'])
    
    # Organized results:
    accuracies = c(accuracies, conf_mtx$overall['Accuracy'])
    
    confusion_matrices[[mdl]] = conf_mtx$table
    
    predictions = cbind(predictions, predict_result)
    colnames(predictions)[dim(predictions)[2]] = paste(mdl, 'predictions', sep='_')
    
    compare_miss = !predict_result==prepared_data$class
    n_misclassifications = rbind(n_misclassifications, t(as.matrix(table(prepared_data$class[compare_miss]))[cell_types, ]))
    
    rownames(n_misclassifications)[dim(n_misclassifications)[1]] = paste(mdl, 'predictions', sep='_')
  }
  names(accuracies) = names(train_results)
  res$accuracies = accuracies
  res$confusion_matrices = confusion_matrices
  res$predictions = predictions
  res$n_misclassifications = n_misclassifications
  return(res)
}

transf_predictions = function(seurat_predict_results){
  res = seurat_predict_results$predictions[,2:dim(seurat_predict_results$predictions)[2]]
  for(iCol in 1:dim(res)[2]){
    res[,iCol] = seurat_predict_results$predictions[, iCol] == seurat_predict_results$predictions$annotated
  }
  return(res)
}

correct_predictions_perCT = function(seurat_predict_results){
  filtered_cells_numbers = c()
  not_filtered_cells_numbers = c()
  transf_preds = transf_predictions(seurat_predict_results)
  for(ct in unique(as.character(seurat_predict_results$predictions$annotated))){
    ct_cells = rownames(seurat_predict_results$predictions)[seurat_predict_results$predictions$annotated==ct]
    filtered_cells_numbers = c(filtered_cells_numbers, sum(rowSums(transf_preds[ct_cells,])>=2))
    not_filtered_cells_numbers = c(not_filtered_cells_numbers, length(ct_cells))
  }
  filtered_cells_numbers = c(filtered_cells_numbers, sum(filtered_cells_numbers))
  not_filtered_cells_numbers = c(not_filtered_cells_numbers, sum(not_filtered_cells_numbers))
  percentages = round(filtered_cells_numbers / not_filtered_cells_numbers * 100, 2)
  res = data.frame(filtered_cells_numbers, not_filtered_cells_numbers, percentages)
  rownames(res) = c(unique(as.character(seurat_predict_results$predictions$annotated)), 'Total')
  colnames(res) = c('Filtered', 'Not Filtered', '%')
  return(res)
}


