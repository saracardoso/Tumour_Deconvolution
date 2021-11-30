code_dir = './code/R'
for(file in list.files(code_dir, full.names=TRUE)) source(file)

# methods = c('AutoGeneS', 'BisqueRNA', 'BSeqSC', 'CIBERSORTx', 'DigitalDLSorter', 'DWLS', 'MOMF', 'MuSiC_woGrouping', 'MuSiC_wGrouping',
#             'Scaden', 'SCDC')
methods = c('AutoGeneS', 'BisqueRNA', 'BSeqSC', 'CIBERSORTx', 'DWLS', 'MOMF', 'MuSiC_woGrouping', 'MuSiC_wGrouping',
            'Scaden', 'SCDC')
cell_types = c('Cancer cells', 'Stromal cells', 'Anti-Inflammatory macro/mono', 'Pro-Inflammatory macro/mono', 'Bcells', 'CD4 Tcells',
               'Regulatory CD4 Tcells', 'CD8 Tcells', 'Proliferative Tcells', 'NK cells', 'Other cells')





# --------------------
# --- TEST RESULTS ---
# --------------------

test_results_dir = './Results/1_Parameter_Optimization/test/proportions'
ground_truth_file = './Data/bulk/CRC/NICs/data/cell_proportions.csv'

methods_colours = c('#ffdb6d', '#c4961a', '#b30000', '#d16103', '#c3d7a4', '#52854c', '#4e84c4','#293352', '#999999', '#cc79a7')
names(methods_colours) = methods



# 1. Get proportions predicted by each method:
predicted_proportions = list()
for(file in list.files(test_results_dir)){
  method_name = gsub('[.]csv', '', file)
  file_path = paste(test_results_dir, file, sep='/')
  predicted_proportions[[method_name]] = read.csv(file_path, row.names=1)
}



# 2. Get ground truth proportions:
ground_truth_proportions = read.csv(ground_truth_file, row.names=1)



# 3. Calculate test rmses:
test_rmses_average_methods = c()
test_rmses_samples_methods = c()
for(method in methods){
  test_rmses_full = rmse_dataset(predicted_proportions[[method]], ground_truth_proportions, cell_types)
  test_rmses_samples_methods = rbind(test_rmses_samples_methods, test_rmses_full$samples)
  test_rmses_average_methods = c(test_rmses_average_methods, test_rmses_full$dataset_average)
}
names(test_rmses_average_methods) = methods
rownames(test_rmses_samples_methods) = methods



# 4. Calculate test correlations:
test_cor_average_methods = c()
test_cor_samples_methods = c()
for(method in methods){
  test_rmses_full = correlation_dataset(predicted_proportions[[method]], ground_truth_proportions, cell_types)
  test_cor_samples_methods = rbind(test_cor_samples_methods, test_rmses_full$samples)
  test_cor_average_methods = c(test_cor_average_methods, test_rmses_full$dataset_average)
}
names(test_cor_average_methods) = methods
rownames(test_cor_samples_methods) = methods



# 5. Scatter plot of rmses vs correlations, averaged through samples:
df_test_scatter_metrics = data.frame(Method=methods,
                                     RMSE=test_rmses_average_methods, Correlation=test_cor_average_methods)
ggplot2::ggplot(df_test_scatter_metrics, ggplot2::aes(RMSE, Correlation, colour=Method)) + 
  ggplot2::geom_point(size=4) + ggplot2::ylim(0,1) + ggplot2::xlim(0,.25) +
  ggplot2::scale_colour_manual(values=c('#ffdb6d', '#c4961a', '#b30000', '#d16103', '#c3d7a4',
                                        '#52854c', '#4e84c4','#293352', '#999999', '#cc79a7'))

# 6. Scatter plot of rmses vs correlations, for each sample:
df_test_scatter_metrics_samples = data.frame(Method=rep(methods, 7),
                                             RMSE=as.numeric(test_rmses_samples_methods), Correlation=as.numeric(test_cor_samples_methods),
                                             CMS=c(rep('CMS2', 10), rep('CMS3', 10), rep('CMS1', 10), rep('Mixed', 10), rep('CMS4', 10),
                                                   rep('Mixed', 10), rep('CMS2', 10)),
                                             Sample=c(rep('NIC5', 10), rep('NIC6', 10), rep('NIC12', 10), rep('NIC16', 10), rep('NIC21', 10),
                                                       rep('NIC27', 10), rep('NIC22nova', 10)))
# 6.1. Coloured by method:
ggplot2::ggplot(df_test_scatter_metrics_samples, ggplot2::aes(RMSE, Correlation, colour=Method)) +
  ggplot2::geom_point(size=2) + ggplot2::xlim(0,.3) +
  ggplot2::scale_colour_manual(values=methods_colours)
# 6.2. Coloured by samples:
ggplot2::ggplot(df_test_scatter_metrics_samples, ggplot2::aes(RMSE, Correlation, colour=Sample)) +
  ggplot2::geom_point(size=2) + ggplot2::xlim(0,.3) +
  ggplot2::scale_colour_manual(values=c('#ffdb6d', '#c4961a', '#b30000', '#d16103', '#c3d7a4',
                                        '#52854c', '#4e84c4','#293352', '#999999', '#cc79a7'))
# 6.3. Coloured by CMS type:
ggplot2::ggplot(df_test_scatter_metrics_samples, ggplot2::aes(RMSE, Correlation, colour=CMS)) +
  ggplot2::geom_point(size=2) + ggplot2::xlim(0,.3) +
  ggplot2::scale_colour_manual(values=c('#b30000', '#d16103', '#52854c', '#4e84c4','#999999', '#cc79a7'))



# 7. Scatter plot of rmses vs correlations, by each cell-type
test_cor_cell_types_methods = c()
test_rmse_cell_types_methods = c()
test_samples = c('NIC5', 'NIC6', 'NIC12', 'NIC16', 'NIC21', 'NIC27', 'NIC22nova')
methods_df = c()
cellTypes_df = c()
for(method in methods){
  for(ct in cell_types){
    methods_df = c(methods_df, method)
    cellTypes_df = c(cellTypes_df, ct)
    test_cor_cell_types_methods = c(test_cor_cell_types_methods, cor(as.numeric(ground_truth_proportions[ct,test_samples]),
                                                                     as.numeric(predicted_proportions[[method]][ct,test_samples]),
                                                                     method = 'pearson'))
    test_rmse_cell_types_methods = c(test_rmse_cell_types_methods,Metrics::rmse(as.numeric(ground_truth_proportions[ct,test_samples]),
                                                                                as.numeric(predicted_proportions[[method]][ct,test_samples])))
  }
}
df_test_scatter_metrics_cellTypes = data.frame(Method=methods_df, RMSE=test_rmse_cell_types_methods, Correlation=test_cor_cell_types_methods,
                                               CellTypes=cellTypes_df)
df_test_scatter_metrics_cellTypes = df_test_scatter_metrics_cellTypes[!is.na(df_test_scatter_metrics_cellTypes$Correlation),]
df_test_scatter_metrics_cellTypes$CellTypes = factor(df_test_scatter_metrics_cellTypes$CellTypes, levels=cell_types)
ggplot2::ggplot(df_test_scatter_metrics_cellTypes, ggplot2::aes(RMSE, Correlation, colour=Method)) + 
  ggplot2::facet_wrap(ggplot2::vars(CellTypes)) +
  ggplot2::geom_point(size=3) + ggplot2::ylim(-.75,1) + ggplot2::xlim(0,.61) +
  ggplot2::scale_colour_manual(values=methods_colours)



# 8. Scatter plots of proportions predicted vs ground truth
plot_estimated_vs_groundTruth_single_method(predicted_proportions$AutoGeneS, ground_truth_proportions, cell_types, cell_type_colors)
plot_estimated_vs_groundTruth_single_method(predicted_proportions$AutoGeneS, ground_truth_proportions, cell_types, cell_type_colors,
                                            annotate=FALSE, xlim=c(0, .2), ylim=c(0, .2))

plot_estimated_vs_groundTruth_single_method(predicted_proportions$BisqueRNA, ground_truth_proportions, cell_types, cell_type_colors)
plot_estimated_vs_groundTruth_single_method(predicted_proportions$BisqueRNA, ground_truth_proportions, cell_types, cell_type_colors,
                                            annotate=FALSE, xlim=c(0, .2), ylim=c(0, .2))

plot_estimated_vs_groundTruth_single_method(predicted_proportions$BSeqSC, ground_truth_proportions, cell_types, cell_type_colors)
plot_estimated_vs_groundTruth_single_method(predicted_proportions$BSeqSC, ground_truth_proportions, cell_types, cell_type_colors,
                                            annotate=FALSE, xlim=c(0, .2), ylim=c(0, .2))

plot_estimated_vs_groundTruth_single_method(predicted_proportions$CIBERSORTx, ground_truth_proportions, cell_types, cell_type_colors)
plot_estimated_vs_groundTruth_single_method(predicted_proportions$CIBERSORTx, ground_truth_proportions, cell_types, cell_type_colors,
                                            annotate=FALSE, xlim=c(0, .2), ylim=c(0, .2))

plot_estimated_vs_groundTruth_single_method(predicted_proportions$DWLS, ground_truth_proportions, cell_types, cell_type_colors)
plot_estimated_vs_groundTruth_single_method(predicted_proportions$DWLS, ground_truth_proportions, cell_types, cell_type_colors,
                                            annotate=FALSE, xlim=c(0, .2), ylim=c(0, .2))

plot_estimated_vs_groundTruth_single_method(predicted_proportions$MOMF, ground_truth_proportions, cell_types, cell_type_colors)
plot_estimated_vs_groundTruth_single_method(predicted_proportions$MOMF, ground_truth_proportions, cell_types, cell_type_colors,
                                            annotate=FALSE, xlim=c(0, .2), ylim=c(0, .2))

plot_estimated_vs_groundTruth_single_method(predicted_proportions$MuSiC_wGrouping, ground_truth_proportions, cell_types, cell_type_colors)
plot_estimated_vs_groundTruth_single_method(predicted_proportions$MuSiC_wGrouping, ground_truth_proportions, cell_types, cell_type_colors,
                                            annotate=FALSE, xlim=c(0, .2), ylim=c(0, .2))

plot_estimated_vs_groundTruth_single_method(predicted_proportions$MuSiC_woGrouping, ground_truth_proportions, cell_types, cell_type_colors)
plot_estimated_vs_groundTruth_single_method(predicted_proportions$MuSiC_woGrouping, ground_truth_proportions, cell_types, cell_type_colors,
                                            annotate=FALSE, xlim=c(0, .2), ylim=c(0, .2))

plot_estimated_vs_groundTruth_single_method(predicted_proportions$Scaden, ground_truth_proportions, cell_types, cell_type_colors)
plot_estimated_vs_groundTruth_single_method(predicted_proportions$Scaden, ground_truth_proportions, cell_types, cell_type_colors,
                                            annotate=FALSE, xlim=c(0, .2), ylim=c(0, .2))

plot_estimated_vs_groundTruth_single_method(predicted_proportions$SCDC, ground_truth_proportions, cell_types, cell_type_colors)
plot_estimated_vs_groundTruth_single_method(predicted_proportions$SCDC, ground_truth_proportions, cell_types, cell_type_colors,
                                            annotate=FALSE, xlim=c(0, .2), ylim=c(0, .2))



# 9. Scatter plots of proportions predicted vs ground truth by method for a cell-type
plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions, ground_truth_proportions, 'Cancer cells', methods_colours,
                                                 xlim=c(0, 1), ylim=c(0, 1))
plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions[c('AutoGeneS', 'CIBERSORTx', 'DWLS', 'MuSiC_woGrouping', 'Scaden')],
                                                 ground_truth_proportions, 'Cancer cells', methods_colours,
                                                 xlim=c(0, 1), ylim=c(0, 1))

plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions, ground_truth_proportions, 'Stromal cells', methods_colours,
                                                 xlim=c(0, 1), ylim=c(0, 1))
plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions[c('BisqueRNA', 'BSeqSC', 'DWLS', 'MuSiC_wGrouping')],
                                                 ground_truth_proportions, 'Stromal cells', methods_colours,
                                                 xlim=c(0, 1), ylim=c(0, 1))

plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions, ground_truth_proportions, 'Anti-Inflammatory macro/mono',
                                                 methods_colours, xlim=c(0, .3), ylim=c(0, .3))
plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions[c('BSeqSC', 'CIBERSORTx', 'MuSiC_wGrouping', 'MuSiC_woGrouping',
                                                                         'Scaden', 'SCDC')],
                                                 ground_truth_proportions, 'Anti-Inflammatory macro/mono', methods_colours,
                                                 xlim=c(0, .1), ylim=c(0, .1))

plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions, ground_truth_proportions, 'Pro-Inflammatory macro/mono',
                                                 methods_colours, xlim=c(0, .2), ylim=c(0, .2))
plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions[c('BisqueRNA', 'DWLS', 'MuSiC_wGrouping', 'MuSiC_woGrouping')],
                                                 ground_truth_proportions, 'Pro-Inflammatory macro/mono', methods_colours,
                                                 xlim=c(0, .2), ylim=c(0, .2))

plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions, ground_truth_proportions, 'Bcells',
                                                 methods_colours, xlim=c(0, .2), ylim=c(0, .2))
plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions[c('CIBERSORTx', 'DWLS', 'MuSiC_woGrouping', 'MuSiC_wGrouping')],
                                                 ground_truth_proportions, 'Bcells', methods_colours,
                                                 xlim=c(0, .15), ylim=c(0, .15))

plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions, ground_truth_proportions, 'CD4 Tcells',
                                                 methods_colours, xlim=c(0, .3), ylim=c(0, .3))
plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions[c('CIBERSORTx', 'AutoGeneS')],
                                                 ground_truth_proportions, 'CD4 Tcells', methods_colours,
                                                 xlim=c(0, .2), ylim=c(0, .2))

plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions, ground_truth_proportions, 'Regulatory CD4 Tcells',
                                                 methods_colours, xlim=c(0, .5), ylim=c(0, .5))
plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions[c('MuSiC_woGrouping', 'AutoGeneS', 'Scaden')],
                                                 ground_truth_proportions, 'Regulatory CD4 Tcells', methods_colours,
                                                 xlim=c(0, .1), ylim=c(0, .1))

plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions, ground_truth_proportions, 'CD8 Tcells',
                                                 methods_colours, xlim=c(0, .5), ylim=c(0, .5))
plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions[c('DWLS', 'AutoGeneS')],
                                                 ground_truth_proportions, 'CD8 Tcells', methods_colours,
                                                 xlim=c(0, .1), ylim=c(0, .1))

plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions, ground_truth_proportions, 'Proliferative Tcells',
                                                 methods_colours, xlim=c(0, .55), ylim=c(0, .55))
plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions[c('DWLS', 'MuSiC_wGrouping', 'Scaden')],
                                                 ground_truth_proportions, 'Proliferative Tcells', methods_colours,
                                                 xlim=c(0, .1), ylim=c(0, .1))

plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions, ground_truth_proportions, 'NK cells',
                                                 methods_colours, xlim=c(0, .3), ylim=c(0, .3))
plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions[c('DWLS', 'BSeqSC', 'Scaden')],
                                                 ground_truth_proportions, 'NK cells', methods_colours,
                                                 xlim=c(0, .1), ylim=c(0, .1))

plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions, ground_truth_proportions, 'NK cells',
                                                 methods_colours, xlim=c(0, .3), ylim=c(0, .3))
plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions[c('DWLS', 'BSeqSC', 'Scaden')],
                                                 ground_truth_proportions, 'NK cells', methods_colours,
                                                 xlim=c(0, .1), ylim=c(0, .1))

plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions, ground_truth_proportions, 'Other cells',
                                                 methods_colours, xlim=c(0, .5), ylim=c(0, .5))
plot_estimated_vs_groundTruth_perMethod_singleCT(predicted_proportions[c('AutoGeneS', 'BisqueRNA', 'BSeqSC', 'CIBERSORTx',
                                                                         'MuSiC_woGrouping', 'SCDC')],
                                                 ground_truth_proportions, 'Other cells', methods_colours,
                                                 xlim=c(0, .3), ylim=c(0, .3))


