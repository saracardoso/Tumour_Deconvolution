code_dir = './code/R'
for(file in list.files(code_dir, full.names=TRUE)) source(file)

methods = c('AutoGeneS_linear', 'AutoGeneS_nnls', 'AutoGeneS_nusvr', 'BisqueRNA', 'BSeqSC', 'CIBERSORTx',
            'DigitalDLSorter', 'DWLS', 'MOMF', 'MuSiC_woGrouping', 'MuSiC_wGrouping', 'Scaden', 'SCDC')
cell_types = c('Cancer cells', 'Stromal cells', 'Macro/mono Lineage', 'Bcells', 'CD4 Tcells',
               'Regulatory CD4 Tcells', 'CD8 Tcells','Proliferative Tcells', 'NK cells')

methods_colours = c('#808000', '#5c4033', '#f2d2bd', '#c4961a', '#b30000', '#d16103', '#ffdb6d',
                    '#c3d7a4','#52854c', '#4e84c4','#293352','#999999', '#cc79a7')
names(methods_colours) = methods

ref_colors = c('#b30000', '#52854c', '#4e84c4')
names(ref_colors) = c('Uncorrected', 'Corrected Before', 'Corrected After')

cell_type_colors = c('#011627', '#8b4513', '#f1558e', '#800080', '#a5be00',
                     '#2ec4b6', '#e71d36', '#3e673c', '#ff9f1c')
names(cell_type_colors) = c('Cancer cells', 'Stromal cells', 'Macro/mono Lineage', 'Bcells', 'CD4 Tcells',
                            'Regulatory CD4 Tcells', 'CD8 Tcells', 'Proliferative Tcells', 'NK cells')







# ---------------------------------------------------
# --- BARPLOT OF GROUNDTRUTH PROPORTIONS / COUNTS ---
# ---------------------------------------------------

# 1. Ground-truth proportions and counts:
ground_truth = read.csv('./Data/bulk/CRC/NICs/data/cell_proportions.csv', row.names=1, check.names = F)[,1:21]
ground_truth_counts = read.csv('./Data/bulk/CRC/NICs/data/cell_counts.csv', row.names=1, check.names=F)[, 1:21]

# 2. Prepare data.frame for plot proportions:
props = c()
cts = c()
samples = c()
for(samp in colnames(ground_truth)){
  props = c(props, ground_truth[, samp])
  cts = c(cts, rownames(ground_truth))
  samples = c(samples, rep(samp, dim(ground_truth)[1]))
}
df_plot = data.frame(props, samples, cts)
cell_type_colors = c(cell_type_colors, 'grey')
names(cell_type_colors)[10] = 'Other cells'
df_plot$cts = factor(df_plot$cts, levels=names(cell_type_colors))

# 3. Plot proportions:
ggplot2::ggplot(df_plot, ggplot2::aes(x=samples, y=props, fill=cts)) +
  ggplot2::geom_col() + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=7)) +
  ggplot2::scale_fill_manual(values=cell_type_colors) + ggplot2::ylab('Proportions')

# 4. Prepare data.frame for plot counts:
counts = c()
cts = c()
samples = c()
for(samp in colnames(ground_truth_counts)){
  counts = c(counts, ground_truth_counts[, samp])
  cts = c(cts, rownames(ground_truth_counts))
  samples = c(samples, rep(samp, dim(ground_truth_counts)[1]))
}
df_plot_counts = data.frame(counts, samples, cts)
cell_type_colors = c(cell_type_colors, 'grey')
names(cell_type_colors)[10] = 'Other cells'
df_plot_counts$cts = factor(df_plot_counts$cts, levels=names(cell_type_colors))

# 5. Plot counts:
ggplot2::ggplot(df_plot_counts, ggplot2::aes(x=samples, y=counts, fill=cts)) +
  ggplot2::geom_col() + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=7)) +
  ggplot2::scale_fill_manual(values=cell_type_colors) + ggplot2::ylab('Cell counts')







# -------------------------------
# --- READ CALCULATED METRICS ---
# -------------------------------

# 0. Metadata:
metadata = data.frame(reference = c('uncorrected', 'correctedBefore', 'correctedAfter', 'uncorrected', 'correctedBefore', 'correctedAfter'),
                      proliferative = c(rep('wProl', 3), rep('woProl', 3), rep('wProl', 3), rep('woProl', 3)),
                      ground_truth = c(rep('gdwOthers', 6), rep('woOthers', 6)),
                      row.names=c('uncorrected_gd', 'correctedBefore_gd', 'correctedAfter_gd',
                                  'uncorrected_woProl_gd', 'correctedBefore_woProl_gd',
                                  'correctedAfter_woProl_gd', 'uncorrected_gdwoOthers',
                                  'correctedBefore_gdwoOthers', 'correctedAfter_gdwoOthers',
                                  'uncorrected_woProl_gdwoOthers', 'correctedBefore_woProl_gdwoOthers',
                                  'correctedAfter_woProl_gdwoOthers'))

# 1. RMSEs:
rmses = as.matrix(read.csv('./Results/metrics/average_rmses.csv', header=TRUE, row.names=1))
rmses_samples = read.csv('./Results/metrics/samples_rmses.csv', header=TRUE, row.names=1)
rmses_cellTypes = read.csv('./Results/metrics/cellTypes_rmses.csv', header=TRUE, row.names=1)

# 2. Correlations:
corrs = as.matrix(read.csv('./Results/metrics/average_correlations.csv', header=TRUE, row.names=1))
corrs_samples = read.csv('./Results/metrics/samples_correlations.csv', header=TRUE, row.names=1)
corrs_cellTypes = read.csv('./Results/metrics/cellTypes_correlations.csv', header=TRUE, row.names=1)





# -------------------------------
# --- PLOT CALCULATED METRICS ---
# -------------------------------


# 1. Heatmaps:
refs = metadata[colnames(rmses),]$reference
refs[grep('uncorrected', refs)] = 'Uncorrected'
refs[grep('correctedBefore', refs)] = 'Corrected Before'
refs[grep('correctedAfter', refs)] = 'Corrected After'

# 1.1. Average RMSEs:
ht1 = ComplexHeatmap::Heatmap(rmses, name='RMSE', heatmap_legend_param = list(direction='horizontal'),
                              show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE,
                        col = circlize::colorRamp2(c(0.1, 0.2, 0.5),
                                                   c('darkgreen', 'white', 'red')),
                        top_annotation = ComplexHeatmap::HeatmapAnnotation(gt=metadata[colnames(rmses),]$ground_truth,
                                                                           prol=metadata[colnames(rmses),]$proliferative,
                                                                           reference=refs,
                                                                           show_annotation_name=FALSE,
                                                                           col=list(reference=ref_colors,
                                                                                    prol=c('wProl'='darkgreen', 'woProl'='darkgrey'),
                                                                                    gt=c('gdwOthers'='gold2', 'woOthers'='darkcyan')),
                                                                           annotation_legend_param=list(gt=list(title='Ground-Truth Proportions',
                                                                                                                at=c('gdwOthers', 'woOthers'),
                                                                                                                labels=c('With Others', 'Without Others')),
                                                                                                        prol=list(title='Proliferative Tcells Predictions',
                                                                                                                  at=c('woProl', 'wProl'),
                                                                                                                  labels=c('Without', 'With')),
                                                                                                        reference=list(title='Reference'))),
                        cell_fun = function(j,i,x,y,width,height,fill){if(!is.na(rmses[i,j])) grid::grid.text(sprintf('%0.3f', rmses[i,j]), x, y,
                                                                                       gp=grid::gpar(fontsize=10))})
ComplexHeatmap::draw(ht1, merge_legend = TRUE, heatmap_legend_side='bottom',
                     annotation_legend_side='bottom')
# 1.1.1. Only with 'with others' and 'with prol'
ht11 = ComplexHeatmap::Heatmap(rmses[,1:3], name='RMSE', heatmap_legend_param=list(direction='horizontal'),
                              show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE,
                              col = circlize::colorRamp2(c(0.1, 0.2, 0.5),
                                                         c('darkgreen', 'white', 'red')),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(reference=refs[1:3],
                                                                                 show_annotation_name=FALSE,
                                                                                 col=list(reference=ref_colors),
                                                                                 annotation_legend_param=list(reference=list(title='Reference'))),
                              cell_fun = function(j,i,x,y,width,height,fill){if(!is.na(rmses[i,j])) grid::grid.text(sprintf('%0.3f', rmses[i,j]), x, y,
                                                                                                                    gp=grid::gpar(fontsize=10))})
ComplexHeatmap::draw(ht11, merge_legend = TRUE, heatmap_legend_side='bottom',
                     annotation_legend_side='bottom')

# 1.2. Average correlations:
ht2 = ComplexHeatmap::Heatmap(corrs, name='p', heatmap_legend_param = list(direction='horizontal'),
                              show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE,
                        col = circlize::colorRamp2(c(-0.5, 0.5, 1),
                                                   c('red', 'white', 'darkgreen')),
                        top_annotation = ComplexHeatmap::HeatmapAnnotation(gt=metadata[colnames(corrs),]$ground_truth,
                                                                           prol=metadata[colnames(corrs),]$proliferative,
                                                                           reference=refs,
                                                                           show_annotation_name=FALSE,
                                                                           col=list(reference=ref_colors,
                                                                                    prol=c('wProl'='darkgreen', 'woProl'='darkgrey'),
                                                                                    gt=c('gdwOthers'='gold2', 'woOthers'='darkcyan')),
                                                                           annotation_legend_param=list(gt=list(title='Ground-Truth Proportions',
                                                                                                                at=c('gdwOthers', 'woOthers'),
                                                                                                                labels=c('With Others', 'Without Others')),
                                                                                                        prol=list(title='Proliferative Tcells Predictions',
                                                                                                                  at=c('woProl', 'wProl'),
                                                                                                                  labels=c('Without', 'With')),
                                                                                                        reference=list(title='Reference'))),
                        cell_fun = function(j,i,x,y,width,height,fill){if(!is.na(corrs[i,j])) grid::grid.text(sprintf('%0.3f', corrs[i,j]), x, y,
                                                                                                              gp=grid::gpar(fontsize=10))})
ComplexHeatmap::draw(ht2, merge_legend = TRUE, heatmap_legend_side='bottom',
                     annotation_legend_side='bottom')
# 1.1.1. Only with 'with others' and 'with prol'
ht21 = ComplexHeatmap::Heatmap(corrs[,1:3], name='p', heatmap_legend_param=list(direction='horizontal'),
                               show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE,
                               col = circlize::colorRamp2(c(-0.5, 0.5, 1),
                                                          c('red', 'white', 'darkgreen')),
                               top_annotation = ComplexHeatmap::HeatmapAnnotation(reference=refs[1:3],
                                                                                  show_annotation_name=FALSE,
                                                                                  col=list(reference=ref_colors),
                                                                                  annotation_legend_param=list(reference=list(title='Reference'))),
                               cell_fun = function(j,i,x,y,width,height,fill){if(!is.na(corrs[i,j])) grid::grid.text(sprintf('%0.3f', corrs[i,j]), x, y,
                                                                                                                     gp=grid::gpar(fontsize=10))})
ComplexHeatmap::draw(ht21, merge_legend = TRUE, heatmap_legend_side='bottom',
                     annotation_legend_side='bottom')


# 2. Scatter plots of samples correlations vs rmses for each rna-correction type, separated by method:
plots_1 = list()
for(cols_type in c('uncorrected_gd', 'correctedBefore_gd', 'correctedAfter_gd')){
  sdf_plta = rmses_samples[rmses_samples$cols==cols_type,]
  colnames(sdf_plta)[1] = 'RMSE'
  sdf_pltb = corrs_samples[rmses_samples$cols==cols_type,]
  colnames(sdf_pltb)[1] = 'Correlation'
  df_plt2 = cbind(sdf_plta, sdf_pltb$Correlation)
  colnames(df_plt2)[5] = 'Correlation'
  averages_df = data.frame(rmse = rmses[, cols_type], correlation =  corrs[, cols_type],
                           models=rownames(rmses))
  plots_1[[cols_type]] = ggplot2::ggplot(df_plt2, ggplot2::aes(RMSE, Correlation, colour=models)) + 
    ggplot2::facet_wrap(ggplot2::vars(models)) +
    ggplot2::geom_point(size=3) + #ggplot2::ylim(-.75,1) + ggplot2::xlim(0,.61) +
    ggplot2::scale_colour_manual(values=methods_colours) +
    ggrepel::geom_text_repel(ggplot2::aes(label=samples), colour='black') +
    ggplot2::geom_vline(xintercept=0.2, linetype='dashed', colour='darkgrey') + 
    ggplot2::geom_hline(yintercept=0.5, linetype='dashed', colour='darkgrey') + 
    ggplot2::geom_point(data=averages_df, mapping=ggplot2::aes(x=rmse, y=correlation),
                        col='black', size=4, shape=8)
}


# 3. Scatter plots of samples correlations vs rmses for each method, separated by rna-correction type:
plots_2 = list()
corrections = c('uncorrected_gd', 'correctedBefore_gd', 'correctedAfter_gd')
for(methods in c("BisqueRNA", "BSeqSC", "CIBERSORTx", "DigitalDLSorter", "DWLS", "MOMF",
                 "MuSiC_wGrouping", "Scaden", "SCDC")){
  sdf_plta = rmses_samples[rmses_samples$models==methods &
                                rmses_samples$cols%in%corrections,]
  colnames(sdf_plta)[1] = 'RMSE'
  sdf_pltb = corrs_samples[rmses_samples$models==methods &
                                       rmses_samples$cols%in%corrections,]
  colnames(sdf_pltb)[1] = 'Correlation'
  df_plt2 = cbind(sdf_plta, sdf_pltb$Correlation)
  colnames(df_plt2)[5] = 'Correlation'
  df_plt2$cols[df_plt2$cols=='uncorrected_gd'] = 'Uncorrected'
  df_plt2$cols[df_plt2$cols=='correctedBefore_gd'] = 'Corrected Before'
  df_plt2$cols[df_plt2$cols=='correctedAfter_gd'] = 'Corrected After'
  df_plt2$cols = factor(df_plt2$cols, levels=names(ref_colors))
  averages_df = data.frame(rmse = rmses[methods, corrections],
                           correlation = corrs[methods, corrections],
                           cols=factor(names(ref_colors), levels=names(ref_colors)))
  plots_2[[methods]] = ggplot2::ggplot(df_plt2, ggplot2::aes(RMSE, Correlation, colour=cols)) + 
    ggplot2::facet_wrap(ggplot2::vars(cols)) +
    ggplot2::geom_point(size=3) + #ggplot2::ylim(-.75,1) + ggplot2::xlim(0,.61) +
    ggplot2::scale_colour_manual(values=ref_colors) +
    ggrepel::geom_text_repel(ggplot2::aes(label=samples), colour='black') +
    ggplot2::geom_vline(xintercept=0.2, linetype='dashed', colour='darkgrey') + 
    ggplot2::geom_hline(yintercept=0.5, linetype='dashed', colour='darkgrey') + 
    ggplot2::geom_point(data=averages_df, mapping=ggplot2::aes(x=rmse, y=correlation),
                        col='black', size=4, shape=8)
}



# 4. Focusing on 'uncorrected' results:

# 4.1. Scatter plots of methods correlations vs rmses, separated by rna-correction type:
df_plt3 = data.frame(Correlation=corrs[,'uncorrected_gd'], RMSE=rmses[,'uncorrected_gd'],
                methods=rownames(corrs))
ggplot2::ggplot(df_plt3, ggplot2::aes(RMSE, Correlation, colour=methods)) +
  ggplot2::geom_point(size=3) + #ggplot2::ylim(-.75,1) + ggplot2::xlim(0,.61) +
  ggplot2::scale_colour_manual(values=methods_colours) +
  ggplot2::geom_vline(xintercept=0.2, linetype='dashed', colour='darkgrey') + 
  ggplot2::geom_hline(yintercept=0.5, linetype='dashed', colour='darkgrey')

# 4.2. Scatter plots of cell-types correlations vs rmses, separated by cell-type and coloured by method:
df_plt3 = cbind(corrs_cellTypes[corrs_cellTypes$cols=='uncorrected_gd', 'values'],
                rmses_cellTypes[corrs_cellTypes$cols=='uncorrected_gd', -4])
colnames(df_plt3)[1:2] = c('Correlation', 'RMSE')
df_plt3$cell_types = factor(df_plt3$cell_types, levels=cell_types)
ggplot2::ggplot(df_plt3, ggplot2::aes(RMSE, Correlation, colour=models)) +
  ggplot2::geom_point(size=3) + #ggplot2::ylim(-.75,1) + ggplot2::xlim(0,.61) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types)) +
  ggplot2::scale_colour_manual(values=methods_colours) +
  ggplot2::geom_vline(xintercept=0.2, linetype='dashed', colour='darkgrey') + 
  ggplot2::geom_hline(yintercept=0.5, linetype='dashed', colour='darkgrey')

# 4.3. Scatter plots of estimated vs ground truth:
# 4.3.1. Read estimated proportions from models' uncorrected, with proliferative Tcells
results_dir = './Results/proportions'
estimated_proportions = list()
res_files = grep('[.]csv', list.files(paste(results_dir, 'uncorrected', sep='/')), value=TRUE)
res_files = res_files[grep('woProl_', res_files, invert=TRUE)]

for(file in res_files){
  res = read.csv(paste(results_dir, 'uncorrected', file, sep='/'), row.names=1)
  model_name =  gsub('[.]csv', '', file)
  estimated_proportions[[model_name]] = res
}
for(model_name in c('AutoGeneS_linear', 'AutoGeneS_nnls', 'AutoGeneS_nusvr', 'MuSiC_woGrouping')){
  file_name = paste(results_dir, '/', model_name, '.csv', sep='')
  res = read.csv(file_name, row.names=1)
  estimated_proportions[[model_name]] = res
}
# Exclude nova samples:
samples_to_keep = c("NIC12", "NIC13", "NIC15", "NIC16", "NIC17", "NIC18", "NIC19", "NIC20", "NIC21", "NIC22", "NIC23",
                    "NIC24", "NIC25", "NIC27", "NIC29", "NIC3", "NIC4", "NIC5", "NIC6", "NIC7", "S52")
for(method in names(estimated_proportions)){
  estimated_proportions[[method]] = estimated_proportions[[method]][,samples_to_keep]
}
# 4.3.2. Read ground-truth proportions:
ground_truth = read.csv('./Data/bulk/CRC/NICs/data/cell_proportions.csv', row.names=1)[, samples_to_keep]
# 4.3.3. Scatter plots of cell-types estimated vs ground truth, separated by method: 
plots_4 = list()
for(model in names(estimated_proportions)){
  plots_4[[model]] = plot_estimated_vs_groundTruth_single_method(estimated_proportions[[model]],
                                                                 ground_truth,
                                                                 cell_types, cell_type_colors,
                                                                 annotate = FALSE)
}
# 4.3.4. Scatter plots of samples estimated vs ground truth for each cell-type, separated by method:
plots_5 = list()
for(ct in cell_types){
  plots_5[[ct]] = plot_estimated_vs_groundTruth_perMethod_singleCT(estimated_proportions,
                                                                   ground_truth, ct, methods_colours)
}
#


# 5. Results combining best methods' results for each cell-type
# 5.1. Combine proportions:
combined_proportions = estimated_proportions$Scaden
combined_proportions['Stromal cells',] = estimated_proportions$BisqueRNA['Stromal cells', colnames(combined_proportions)]
combined_proportions['Macro/mono Lineage',] = estimated_proportions$BSeqSC['Macro/mono Lineage', colnames(combined_proportions)]
combined_proportions['CD4 Tcells',] = estimated_proportions$CIBERSORTx['CD4 Tcells', colnames(combined_proportions)]
combined_proportions['CD8 Tcells',] = estimated_proportions$DigitalDLSorter['CD8 Tcells', colnames(combined_proportions)]
combined_proportions['Proliferative Tcells',] = estimated_proportions$MuSiC_wGrouping['Proliferative Tcells', colnames(combined_proportions)]
# 5.2. Proportions normalised for sum to 1:
combined_proportions_norm = combined_proportions
for(samp in colnames(combined_proportions)){
  sum_samp = sum(combined_proportions[,samp])
  combined_proportions_norm[,samp] = combined_proportions[,samp] / sum_samp
}
# 5.3. correlation and rmse:
rmse_comb = rmse_dataset(combined_proportions, ground_truth, rownames(combined_proportions))
corr_comb = correlation_dataset(combined_proportions, ground_truth, rownames(combined_proportions))
rmse_comb_norm = rmse_dataset(combined_proportions_norm, ground_truth, rownames(combined_proportions))
corr_comb_norm = correlation_dataset(combined_proportions_norm, ground_truth, rownames(combined_proportions))
# 5.4. Estimated vs ground-truth proportions
plot_estimated_vs_groundTruth_single_method(combined_proportions, ground_truth, cell_types, cell_type_colors,
                                            annotate = TRUE)
plot_estimated_vs_groundTruth_single_method(combined_proportions_norm, ground_truth, cell_types, cell_type_colors,
                                            annotate = TRUE)
# 5.5. Correlation vs rmse:
corr = c(corrs[,'uncorrected_gd'], corr_comb$dataset_average, corr_comb_norm$dataset_average)
names(corr) = c(rownames(corrs), 'Combined', 'Combined_norm')
rmse = c(rmses[,'uncorrected_gd'], rmse_comb$dataset_average, rmse_comb_norm$dataset_average)
methods_colours_comb = c(methods_colours, '#000000', '#b2beb5')
names(methods_colours_comb) = c(names(methods_colours), 'Combined', 'Combined_norm')
df_plt_comb = data.frame(Correlation=corr, RMSE=rmse, Methods=names(corr),
                         Type=c(rep('Original', 13), 'Combined', 'Combined'))
df_plt_comb$Type = factor(df_plt_comb$Type, levels=c('Original', 'Combined'))
ggplot2::ggplot(df_plt_comb, ggplot2::aes(RMSE, Correlation, colour=Methods, shape=Type)) +
  ggplot2::geom_point(size=4) + #ggplot2::ylim(-.75,1) + ggplot2::xlim(0,.61) +
  ggplot2::scale_colour_manual(values=methods_colours_comb) +
  ggplot2::scale_shape_manual(values=c('Original'=19, 'Combined'=18)) +
  ggplot2::geom_vline(xintercept=0.2, linetype='dashed', colour='darkgrey') + 
  ggplot2::geom_hline(yintercept=0.5, linetype='dashed', colour='darkgrey')
# 5.6. correlation vs rmse per sample combined vs not normalised vs scaden:
scaden_corrs = corrs_samples[corrs_samples$models=='Scaden' & corrs_samples$cols=='uncorrected_gd', 'values']
names(scaden_corrs) = corrs_samples[corrs_samples$models=='Scaden' & corrs_samples$cols=='uncorrected_gd', 'samples']
scaden_rmses = rmses_samples[rmses_samples$models=='Scaden' & rmses_samples$cols=='uncorrected_gd', 'values']
names(scaden_rmses) = rmses_samples[rmses_samples$models=='Scaden' & rmses_samples$cols=='uncorrected_gd', 'samples']

df_plt_comb_samp = data.frame(Correlation=c(corr_comb$samples, corr_comb_norm$samples, scaden_corrs[names(corr_comb_norm$samples)]),
                              RMSE=c(rmse_comb$samples, rmse_comb_norm$samples, scaden_rmses[names(corr_comb_norm$samples)]),
                              Method=c(rep('Combined', length(corr_comb$samples)), rep('Norm Combined', length(corr_comb$samples)),
                                       rep('Scaden', length(corr_comb$samples))),
                              samples = c(names(corr_comb$samples), names(rmse_comb$samples), names(rmse_comb$samples)))
ggplot2::ggplot(df_plt_comb_samp, ggplot2::aes(RMSE, Correlation, colour=Method)) +
  ggplot2::geom_point(size=3) + ggplot2::facet_wrap(ggplot2::vars(samples), ncol=6) +
  #ggrepel::geom_text_repel(ggplot2::aes(label=Method), colour='black') +
  ggplot2::scale_colour_manual(values=c('Scaden'='#999999', 'Combined'='#52854c', 'Norm Combined'='#4e84c4')) +
  ggplot2::geom_vline(xintercept=0.2, linetype='dashed', colour='darkgrey') + 
  ggplot2::geom_hline(yintercept=0.5, linetype='dashed', colour='darkgrey')
  
  
# 6. Check RNA content bias correction per cell-type
x = expand.grid(c('uncorrected', 'after', 'before'), cell_types)
x_names = paste(x$Var2, x$Var1, sep='_')
rmses_cts = matrix(data=NA, nrow=27, ncol=9,
                   dimnames=list(x_names,
                                 c('BisqueRNA', 'BSeqSC', 'CIBERSORTx', 'DigitalDLSorter', 'DWLS', 'MOMF', 'MuSiC_wGrouping',
                                   'Scaden', 'SCDC')))
corrs_cts = matrix(data=NA, nrow=27, ncol=9,
                   dimnames=list(x_names,
                                 c('BisqueRNA', 'BSeqSC', 'CIBERSORTx', 'DigitalDLSorter', 'DWLS', 'MOMF', 'MuSiC_wGrouping',
                                   'Scaden', 'SCDC')))

corrections = c('uncorrected_gd', 'correctedBefore_gd', 'correctedAfter_gd')
for(ct in cell_types){
  for(correction in corrections){
   if(correction == 'uncorrected_gd') corr2 = 'uncorrected'
   else if(correction == 'correctedBefore_gd') corr2 = 'before'
   else corr2 = 'after'
    rmses_ctcor = rmses_cellTypes[rmses_cellTypes$cell_types==ct &
                               rmses_cellTypes$cols==correction,]
    rownames(rmses_ctcor) = rmses_ctcor$models
    rmses_cts[paste(ct, corr2, sep='_'),] = rmses_ctcor[colnames(rmses_cts), 'values']
    
    corrs_ctcor = corrs_cellTypes[corrs_cellTypes$cell_types==ct &
                                    corrs_cellTypes$cols==correction,]
    rownames(corrs_ctcor) = corrs_ctcor$models
    corrs_cts[paste(ct, corr2, sep='_'),] = corrs_ctcor[colnames(corrs_cts), 'values']
  }
}
# 6.1. Heatmap of RMSEs
ht6 = ComplexHeatmap::Heatmap(rmses_cts, name='p', heatmap_legend_param = list(direction='horizontal'),
                              show_column_names = FALSE, cluster_columns = FALSE, show_row_names = FALSE, cluster_rows = FALSE,
                              col = circlize::colorRamp2(c(0.1, 0.2, 0.5),
                                                         c('darkgreen', 'white', 'red')),
                              row_split = x$Var2, row_title=NULL, 
                              left_annotation = ComplexHeatmap::rowAnnotation(cell_type=x$Var2,correction=x$Var1,
                                                                               show_annotation_name=FALSE,
                                                                               col=list(correction=c('uncorrected'="#b30000", 
                                                                                                     'before'="#52854c",
                                                                                                     'after'= "#4e84c4"),
                                                                                        cell_type=cell_type_colors),
                                                                               annotation_legend_param=list(cell_type=list(title='Cell-type', ncol=1),
                                                                                                            correction=list(title='Correction', ncol=2))),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(methods=colnames(rmses_cts), show_annotation_name=FALSE,
                                                                                 col=list(methods=methods_colours[colnames(rmses_cts)]),
                                                                                 annotation_legend_param=list(methods=list(title='Methods', ncol=1))),
                              cell_fun = function(j,i,x,y,width,height,fill){if(!is.na(rmses_cts[i,j])) grid::grid.text(sprintf('%0.3f', rmses_cts[i,j]), x, y,
                                                                                                                    gp=grid::gpar(fontsize=10))})
ComplexHeatmap::draw(ht6, merge_legend = TRUE, heatmap_legend_side='right',
                     annotation_legend_side='right', 
                     column_title='Methods', row_title='Correction Type per Cell-type', row_title_side='left')
# 6.2. Heatmap of correlations
ht62 = ComplexHeatmap::Heatmap(corrs_cts, name='p', heatmap_legend_param = list(direction='horizontal'),
                              show_column_names = FALSE, cluster_columns = FALSE, show_row_names = FALSE, cluster_rows = FALSE,
                              col = circlize::colorRamp2(c(-0.5, 0.5, 1),
                                                         c('red', 'white', 'darkgreen')),
                              row_split = x$Var2, row_title=NULL,
                              left_annotation = ComplexHeatmap::rowAnnotation(correction=x$Var1,cell_type=x$Var2,
                                                                               show_annotation_name=FALSE,
                                                                               col=list(correction=c('uncorrected'="#b30000", 
                                                                                                     'before'="#52854c",
                                                                                                     'after'= "#4e84c4"),
                                                                                        cell_type=cell_type_colors),
                                                                               annotation_legend_param=list(cell_type=list(title='Cell-type', ncol=1),
                                                                                                            correction=list(title='Correction', ncol=2))),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(methods=colnames(corrs_cts), show_annotation_name=FALSE,
                                                                                 col=list(methods=methods_colours[colnames(corrs_cts)]),
                                                                                 annotation_legend_param=list(methods=list(title='Methods', ncol=1))),
                              cell_fun = function(j,i,x,y,width,height,fill){if(!is.na(corrs_cts[i,j])) grid::grid.text(sprintf('%0.3f', corrs_cts[i,j]), x, y,
                                                                                                                        gp=grid::gpar(fontsize=10))})
ComplexHeatmap::draw(ht62, merge_legend = TRUE, heatmap_legend_side='right',
                     annotation_legend_side='right',
                     column_title='Methods', row_title='Correction Type per Cell-type', row_title_side='left')
# 6.3. Scatter plots of estimations vs ground-truth
# 6.3.1. Read results:
results_dir = './Results/proportions'
samples_to_keep = c("NIC12", "NIC13", "NIC15", "NIC16", "NIC17", "NIC18", "NIC19", "NIC20", "NIC21", "NIC22", "NIC23",
                    "NIC24", "NIC25", "NIC27", "NIC29", "NIC3", "NIC4", "NIC5", "NIC6", "NIC7", "S52")
estimated_proportions_bias = list()
for(method in c("BisqueRNA", "BSeqSC", "CIBERSORTx", "DigitalDLSorter", "DWLS", "MOMF", "MuSiC_wGrouping", "Scaden", "SCDC"))
  estimated_proportions_bias[[method]] = list()
for(correction in c('uncorrected', 'corrected_before', 'corrected_after')){
  res_files = grep('[.]csv', list.files(paste(results_dir, correction, sep='/')), value=TRUE)
  res_files = res_files[grep('woProl_', res_files, invert=TRUE)]
  
  for(file in res_files){
    res = read.csv(paste(results_dir, correction, file, sep='/'), row.names=1)
    model_name =  gsub('[.]csv', '', file)
    estimated_proportions_bias[[model_name]][[correction]] = res[, samples_to_keep]
  }
}
# 6.3.1. Plot for each method:
# Cancer cells
props = list()
props[['Scaden Uncorrected']] = estimated_proportions[['Scaden']]
props[['Scaden Corrected Before']] = estimated_proportions_bias$Scaden$corrected_before
plot_estimated_vs_groundTruth_perMethod_singleCT(props, ground_truth, 'Cancer cells', methods_colours, scales='fixed')
# Stromal cells
props = list()
props[['BisqueRNA Uncorrected']] = estimated_proportions[['BisqueRNA']]
props[['BisqueRNA Corrected After']] = estimated_proportions_bias$BisqueRNA$corrected_after
props[['MOMF Corrected Before']] = estimated_proportions_bias$MOMF$corrected_before
props[['MuSiC_wGrouping Corrected After']] = estimated_proportions_bias$MuSiC_wGrouping$corrected_after
plot_estimated_vs_groundTruth_perMethod_singleCT(props, ground_truth, 'Stromal cells', methods_colours, scales='fixed')
# Macro/mono Lineage
props = list()
props[['BSeqSC Uncorrected']] = estimated_proportions[['BSeqSC']]
props[['BisqueRNA Corrected After']] = estimated_proportions_bias$BisqueRNA$corrected_after
props[['BSeqSC Corrected Before']] = estimated_proportions_bias$BSeqSC$corrected_before
plot_estimated_vs_groundTruth_perMethod_singleCT(props, ground_truth, 'Macro/mono Lineage', methods_colours, scales='fixed')
# ---- Bcells
props = list()
props[['Scaden Uncorrected']] = estimated_proportions[['Scaden']]
props[['DigitalDLSorter Corrected Before']] = estimated_proportions_bias$DigitalDLSorter$corrected_before
plot_estimated_vs_groundTruth_perMethod_singleCT(props, ground_truth, 'Bcells', methods_colours, scales='fixed')
# ---- CD4 T-cells
props = list()
props[['CIBERSORTx Uncorrected']] = estimated_proportions[['CIBERSORTx']]
props[['BSeqSC Corrected Before']] = estimated_proportions_bias$BSeqSC$corrected_before
props[['MuSiC_wGrouping Corrected After']] = estimated_proportions_bias$MuSiC_wGrouping$corrected_after
plot_estimated_vs_groundTruth_perMethod_singleCT(props, ground_truth, 'CD4 Tcells', methods_colours, scales='fixed')
# ---- Regulatory CD4 T-cells
props = list()
props[['Scaden Uncorrected']] = estimated_proportions[['Scaden']]
props[['BSeqSC Corrected Before']] = estimated_proportions_bias$BSeqSC$corrected_before
plot_estimated_vs_groundTruth_perMethod_singleCT(props, ground_truth, 'Regulatory CD4 Tcells', methods_colours, scales='fixed')
# ---- CD8 T-cells
props = list()
props[['DigitalDLSorter Uncorrected']] = estimated_proportions[['DigitalDLSorter']]
props[['MOMF Corrected Before']] = estimated_proportions_bias$MOMF$corrected_before
plot_estimated_vs_groundTruth_perMethod_singleCT(props, ground_truth, 'CD8 Tcells', methods_colours, scales='fixed')
# ---- Proliferative T-cells
props = list()
props[['MuSiC_wGrouping Uncorrected']] = estimated_proportions[['MuSiC_wGrouping']]
props[['BisqueRNA Corrected After']] = estimated_proportions_bias$BisqueRNA$corrected_after
props[['DigitalDLSorter Corrected Before']] = estimated_proportions_bias$BisqueRNA$corrected_before
plot_estimated_vs_groundTruth_perMethod_singleCT(props, ground_truth, 'Proliferative Tcells', methods_colours, scales='fixed')
# ---- NK cells
props = list()
props[['Scaden Uncorrected']] = estimated_proportions[['Scaden']]
props[['DigitalDLSorter Corrected Before']] = estimated_proportions_bias$DigitalDLSorter$corrected_before
plot_estimated_vs_groundTruth_perMethod_singleCT(props, ground_truth, 'NK cells', methods_colours, scales='fixed')





# ---------------------------------------------------------------
# --- GET INFO ABOUT CELL COUNTS, READ COUNTS AND PROPORTIONS ---
# ---------------------------------------------------------------

# 1. Get cell counts info:
cell_counts = read.csv('./Data/bulk/CRC/NICs/data/cell_counts.csv', row.names=1)[,samples_to_keep]
total_cell_counts = colSums(cell_counts)

# 2. Get cell counts info:
bulk_data = read.csv('./Data/bulk/CRC/NICs/data/data.csv', row.names=1)[,samples_to_keep]
total_read_counts = colSums(bulk_data)

# 3. Get cell props info:
cell_props = read.csv('./Data/bulk/CRC/NICs/data/cell_proportions.csv', row.names=1)[,samples_to_keep]
total_cell_props = colSums(cell_props)

# 4. Cell VS read counts:
df_cell_read = data.frame(cell_counts=total_cell_counts, read_counts=total_read_counts[names(total_cell_counts)],
                          samples=names(total_cell_counts))
p_corr = paste('p =', round(cor(df_cell_read$cell_counts, df_cell_read$read_counts), 3))
ggplot2::ggplot(df_cell_read, ggplot2::aes(cell_counts, read_counts)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method='lm', colour='#80b3ff') +
  ggplot2::geom_text(data=data.frame(x=12000, y=10000000, lab=p_corr),
                     ggplot2::aes(x, y, label=lab), parse=FALSE, color='black', size=3) +
  ggrepel::geom_text_repel(ggplot2::aes(label=samples), colour='black') +
  ggplot2::theme(legend.position = 'none') + ggplot2::ylab('Total Read Counts') + ggplot2::xlab('Total Cell counts')





# -----------------------------------------------------
# --- COMPARE PREDICTION ABILITY TO NUMBER OF CELLS ---
# -----------------------------------------------------

# 1. Base dataframes for plots:
# 1.1. RMSEs
df_rmse = rmses_samples[rmses_samples$cols=='uncorrected_gd',]
df_rmse$total_counts = 0
df_rmse$cancer_counts = 0
df_rmse$stromal_counts = 0
df_rmse$immune_counts = 0
df_rmse$others_counts = 0
for(samp in names(total_cell_counts)){
  df_rmse[df_rmse$samples==samp, 'total_counts'] = total_cell_counts[samp]
  df_rmse[df_rmse$samples==samp, 'cancer_counts'] = sum(cell_counts[1, samp])
  df_rmse[df_rmse$samples==samp, 'stromal_counts'] = sum(cell_counts[2, samp])
  df_rmse[df_rmse$samples==samp, 'immune_counts'] = sum(cell_counts[-c(1,2,10), samp])
  df_rmse[df_rmse$samples==samp, 'others_counts'] = sum(cell_counts[10, samp])
}

# 1.2. Correlations
df_corr = corrs_samples[corrs_samples$cols=='uncorrected_gd',]
df_corr$total_counts = 0
df_corr$cancer_counts = 0
df_corr$stromal_counts = 0
df_corr$immune_counts = 0
df_corr$others_counts = 0
for(samp in names(total_cell_counts)){
  df_corr[df_corr$samples==samp, 'total_counts'] = total_cell_counts[samp]
  df_corr[df_corr$samples==samp, 'cancer_counts'] = sum(cell_counts[1, samp])
  df_corr[df_corr$samples==samp, 'stromal_counts'] = sum(cell_counts[2, samp])
  df_corr[df_corr$samples==samp, 'immune_counts'] = sum(cell_counts[-c(1,2,10), samp])
  df_corr[df_corr$samples==samp, 'others_counts'] = sum(cell_counts[10, samp])
}

# 2. Add labels:
df_rmse$label = df_rmse$samples
df_corr$label = df_corr$samples


# 2. Plot how samples rmses correlate to counts, separated by method: (focus on the 5 best models)
#models = c('CIBERSORTx', 'DigitalDLSorter', 'Scaden', 'AutoGeneS_nusvr', 'DWLS')#unique(df_rmse$models))
rmses_plots = list()
for(comparison_type in c('total_counts', 'cancer_counts', 'stromal_counts',
                         'immune_counts', 'others_counts')){
  fited_res = c()
  min_rmse = c()
  min_counts = c()
  for(model in methods){
    fited_res = c(fited_res, paste('p =',
                                   round(cor(df_rmse$values[df_rmse$models==model],
                                       df_rmse[df_rmse$models==model, comparison_type]), 3)))
    min_rmse = c(min_rmse, min(df_rmse[df_rmse$models==model,'values']))
    min_counts = c(min_counts,
                   max(df_rmse[df_rmse$models==model, comparison_type]) -
                     (max(df_rmse[df_rmse$models==model, comparison_type])/10))
  }
  names(fited_res) = methods
  
  rmses_plots[[comparison_type]] = ggplot2::ggplot(df_rmse[df_rmse$models%in%methods,],
                                    ggplot2::aes_string(comparison_type, 'values', colour='models')) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(ggplot2::vars(models), scales='free') +
    #ggrepel::geom_text_repel(ggplot2::aes(label=label), colour='black') +
    ggplot2::scale_colour_manual(values=methods_colours) +
    ggplot2::geom_hline(yintercept=0.2, linetype='dashed', colour='darkgrey') +
    ggplot2::geom_smooth(method='lm', colour='#80b3ff') +
    ggplot2::geom_text(data=data.frame(x=min_counts, y=min_rmse, lab=fited_res, models=names(fited_res)),
                       ggplot2::aes(x, y, label=lab), parse=FALSE, color='black', size=3) +
    ggplot2::theme(legend.position = 'none') + ggplot2::ylab('RMSE') + ggplot2::xlab('Cell counts')
}
#
rmses_plots$total_counts
rmses_plots$cancer_counts

# 3. Plot how samples correlations correlate to counts, separated by method:
#models = c('CIBERSORTx', 'DigitalDLSorter', 'Scaden', 'AutoGeneS_nusvr', 'DWLS')#unique(df_rmse$models))
corrs_plots = list()
for(comparison_type in c('total_counts', 'cancer_counts', 'stromal_counts',
                         'immune_counts', 'others_counts')){
  fited_res = c()
  max_corr = c()
  min_counts = c()
  for(model in methods){
    fited_res = c(fited_res, paste('p =',
                                   round(cor(df_corr$values[df_corr$models==model],
                                             df_corr[df_corr$models==model, comparison_type]), 3)))
    max_corr = c(max_corr, max(df_corr[df_corr$models==model,'values']))
    min_counts = c(min_counts,
                   max(df_corr[df_corr$models==model, comparison_type]) -
                     (max(df_corr[df_corr$models==model, comparison_type])/10))
  }
  names(fited_res) = methods
  
  corrs_plots[[comparison_type]] = ggplot2::ggplot(df_corr[df_corr$models%in%methods,],
                                                   ggplot2::aes_string(comparison_type, 'values',
                                                                       colour='models')) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(ggplot2::vars(models), scales='free') +
    ggrepel::geom_text_repel(ggplot2::aes(label=label), colour='black') +
    ggplot2::scale_colour_manual(values=methods_colours) +
    ggplot2::geom_hline(yintercept=0.5, linetype='dashed', colour='darkgrey') +
    ggplot2::geom_smooth(method='lm', colour='#80b3ff') +
    ggplot2::geom_text(data=data.frame(x=min_counts, y=max_corr, lab=fited_res, models=names(fited_res)),
                       ggplot2::aes(x, y, label=lab), parse=FALSE, color='black', size=3) +
    ggplot2::theme(legend.position = 'none') + ggplot2::ylab('Correlation') + ggplot2::xlab('Cell counts')
}
#





# -------------------------------------------------------
# --- COMPARE PREDICTION ABILITY TO TOTAL READ COUNTS ---
# -------------------------------------------------------

# 1. Base dataframes for plots:
# 1.1. RMSEs
df_rmse_rc = rmses_samples[rmses_samples$cols=='uncorrected_gd',]
df_rmse_rc$read_counts = total_read_counts

# 1.2. Correlations
df_corr_rc = corrs_samples[corrs_samples$cols=='uncorrected_gd',]
df_corr_rc$read_counts = total_read_counts


# 2. Add labels:
df_rmse_rc$label = df_rmse_rc$samples
df_corr_rc$label = df_corr_rc$samples


# 3. Plot how samples rmses correlate to counts, separated by method: (focus on the 5 best models)
#models = c('CIBERSORTx', 'DigitalDLSorter', 'Scaden', 'AutoGeneS_nusvr', 'DWLS')#unique(df_rmse_rc$models))
fited_res = c()
min_rmse = c()
min_counts = c()
for(model in methods){
  fited_res = c(fited_res, paste('p =',
                                 round(cor(df_rmse_rc$values[df_rmse_rc$models==model],
                                           df_rmse_rc[df_rmse_rc$models==model, 'read_counts']), 3)))
  min_rmse = c(min_rmse, min(df_rmse_rc[df_rmse_rc$models==model,'values']))
  min_counts = c(min_counts,
                 max(df_rmse_rc[df_rmse_rc$models==model, 'read_counts']) -
                   (max(df_rmse_rc[df_rmse_rc$models==model, 'read_counts'])/10))
}
names(fited_res) = methods

ggplot2::ggplot(df_rmse_rc[df_rmse_rc$models%in%methods,],
                ggplot2::aes_string('read_counts', 'values', colour='models')) +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(ggplot2::vars(models), scales='free') +
  #ggrepel::geom_text_repel(ggplot2::aes(label=label), colour='black') +
  ggplot2::scale_colour_manual(values=methods_colours) +
  ggplot2::geom_hline(yintercept=0.2, linetype='dashed', colour='darkgrey') +
  ggplot2::geom_smooth(method='lm', colour='#80b3ff') +
  ggplot2::geom_text(data=data.frame(x=min_counts, y=min_rmse, lab=fited_res, models=names(fited_res)),
                     ggplot2::aes(x, y, label=lab), parse=FALSE, color='black', size=3) +
  ggplot2::theme(legend.position = 'none') + ggplot2::ylab('RMSE') + ggplot2::xlab('Read counts')
#


# 4. Plot how samples correlations correlate to counts, separated by method: (focus on the 5 best models)
models = c('CIBERSORTx', 'DigitalDLSorter', 'Scaden', 'AutoGeneS_nusvr', 'DWLS')#unique(df_corr_rc$models))
fited_res = c()
min_rmse = c()
min_counts = c()
for(model in models){
  fited_res = c(fited_res, paste('p =',
                                 round(cor(df_corr_rc$values[df_corr_rc$models==model],
                                           df_corr_rc[df_corr_rc$models==model, 'read_counts']), 3)))
  min_rmse = c(min_rmse, min(df_corr_rc[df_corr_rc$models==model,'values']))
  min_counts = c(min_counts,
                 max(df_corr_rc[df_corr_rc$models==model, 'read_counts']) -
                   (max(df_corr_rc[df_corr_rc$models==model, 'read_counts'])/10))
}
names(fited_res) = models

ggplot2::ggplot(df_corr_rc[df_corr_rc$models%in%models,],
                ggplot2::aes_string('read_counts', 'values', colour='models')) +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(ggplot2::vars(models), scales='free') +
  ggrepel::geom_text_repel(ggplot2::aes(label=label), colour='black') +
  ggplot2::scale_colour_manual(values=methods_colours) +
  ggplot2::geom_hline(yintercept=0.5, linetype='dashed', colour='darkgrey') +
  ggplot2::geom_smooth(method='lm', colour='#80b3ff') +
  ggplot2::geom_text(data=data.frame(x=min_counts, y=min_rmse, lab=fited_res, models=names(fited_res)),
                     ggplot2::aes(x, y, label=lab), parse=FALSE, color='black', size=3) +
  ggplot2::theme(legend.position = 'none') + ggplot2::ylab('Correlation') + ggplot2::xlab('Cell counts')
#





# ----------------------------------------------------------
# --- COMPARE PREDICTION ABILITY TO PROPORTIONS OF CELLS ---
# ----------------------------------------------------------

# 1. Base dataframes for plots:
# 1.1. RMSEs
df_rmse_p = rmses_samples[rmses_samples$cols=='uncorrected_gd',]
df_rmse_p$cancer_props = 0
df_rmse_p$stromal_props = 0
df_rmse_p$immune_props = 0
df_rmse_p$others_props = 0
for(samp in names(total_cell_props)){
  df_rmse_p[df_rmse_p$samples==samp, 'cancer_props'] = sum(cell_props[1, samp])
  df_rmse_p[df_rmse_p$samples==samp, 'stromal_props'] = sum(cell_props[2, samp])
  df_rmse_p[df_rmse_p$samples==samp, 'immune_props'] = sum(cell_props[-c(1,2,10), samp])
  df_rmse_p[df_rmse_p$samples==samp, 'others_props'] = sum(cell_props[10, samp])
}

# 1.2. Correlations
df_corr_p = corrs_samples[corrs_samples$cols=='uncorrected_gd',]
df_corr_p$total_props = 0
df_corr_p$cancer_props = 0
df_corr_p$stromal_props = 0
df_corr_p$immune_props = 0
df_corr_p$others_props = 0
for(samp in names(total_cell_counts)){
  df_corr_p[df_corr_p$samples==samp, 'cancer_props'] = sum(cell_props[1, samp])
  df_corr_p[df_corr_p$samples==samp, 'stromal_props'] = sum(cell_props[2, samp])
  df_corr_p[df_corr_p$samples==samp, 'immune_props'] = sum(cell_props[-c(1,2,10), samp])
  df_corr_p[df_corr_p$samples==samp, 'others_props'] = sum(cell_props[10, samp])
}

# 2. Add labels:
df_rmse_p$label = df_rmse_p$samples
df_corr_p$label = df_rmse_p$samples


# 2. Plot how samples rmses correlate to counts, separated by method: (focus on the 5 best models)
#models = c('CIBERSORTx', 'DigitalDLSorter', 'Scaden', 'AutoGeneS_nusvr', 'DWLS') #unique(df_rmse$models))
rmses_plots = list()
for(comparison_type in c('cancer_props', 'stromal_props', 'immune_props', 'others_props')){
  fited_res = c()
  min_rmse = c()
  min_counts = c()
  for(model in methods){
    fited_res = c(fited_res, paste('p =',
                                   round(cor(df_rmse_p$values[df_rmse_p$models==model],
                                             df_rmse_p[df_rmse_p$models==model, comparison_type]), 3)))
    min_rmse = c(min_rmse, min(df_rmse_p[df_rmse_p$models==model,'values']))
    min_counts = c(min_counts,
                   max(df_rmse_p[df_rmse_p$models==model, comparison_type]) -
                     (max(df_rmse_p[df_rmse_p$models==model, comparison_type])/10))
  }
  names(fited_res) = methods
  
  rmses_plots[[comparison_type]] = ggplot2::ggplot(df_rmse_p[df_rmse_p$models%in%methods,],
                                                   ggplot2::aes_string(comparison_type, 'values', colour='models')) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(ggplot2::vars(models), scales='free') +
    #ggrepel::geom_text_repel(ggplot2::aes(label=label), colour='black') +
    ggplot2::scale_colour_manual(values=methods_colours) +
    ggplot2::geom_hline(yintercept=0.2, linetype='dashed', colour='darkgrey') +
    ggplot2::geom_smooth(method='lm', colour='#80b3ff') +
    ggplot2::geom_text(data=data.frame(x=min_counts, y=min_rmse, lab=fited_res, models=names(fited_res)),
                       ggplot2::aes(x, y, label=lab), parse=FALSE, color='black', size=3) +
    ggplot2::theme(legend.position = 'none') + ggplot2::ylab('RMSE') + ggplot2::xlab('Cell counts')
}
#
rmses_plots$cancer_props
rmses_plots$stromal_props
rmses_plots$others_props
rmses_plots$immune_props

# 3. Plot how samples correlations correlate to counts, separated by method:
#models = c('CIBERSORTx', 'DigitalDLSorter', 'Scaden', 'AutoGeneS_nusvr', 'DWLS')#unique(df_corr_p$models))
corrs_plots = list()
for(comparison_type in c('cancer_props', 'stromal_props', 'immune_props', 'others_props')){
  fited_res = c()
  max_corr = c()
  min_counts = c()
  for(model in models){
    fited_res = c(fited_res, paste('p =',
                                   round(cor(df_corr_p$values[df_corr_p$models==model],
                                             df_corr_p[df_corr_p$models==model, comparison_type]), 3)))
    max_corr = c(max_corr, max(df_corr_p[df_corr_p$models==model,'values']))
    min_counts = c(min_counts,
                   max(df_corr_p[df_corr_p$models==model, comparison_type]) -
                     (max(df_corr_p[df_corr_p$models==model, comparison_type])/10))
  }
  names(fited_res) = models
  
  corrs_plots[[comparison_type]] = ggplot2::ggplot(df_corr_p[df_corr_p$models%in%models,],
                                                   ggplot2::aes_string(comparison_type, 'values',
                                                                       colour='models')) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(ggplot2::vars(models), scales='free') +
    ggrepel::geom_text_repel(ggplot2::aes(label=label), colour='black') +
    ggplot2::scale_colour_manual(values=methods_colours) +
    ggplot2::geom_hline(yintercept=0.5, linetype='dashed', colour='darkgrey') +
    ggplot2::geom_smooth(method='lm', colour='#80b3ff') +
    ggplot2::geom_text(data=data.frame(x=min_counts, y=max_corr, lab=fited_res, models=names(fited_res)),
                       ggplot2::aes(x, y, label=lab), parse=FALSE, color='black', size=3) +
    ggplot2::theme(legend.position = 'none') + ggplot2::ylab('Correlation') + ggplot2::xlab('Cell counts')
}
#

