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
rmses = as.matrix(read.csv('./Results/2_Default_Parameters/metrics/average_rmses.csv', header=TRUE, row.names=1))
rmses_samples = read.csv('./Results/2_Default_Parameters/metrics/samples_rmses.csv', header=TRUE, row.names=1)
rmses_cellTypes = read.csv('./Results/2_Default_Parameters/metrics/cellTypes_rmses.csv', header=TRUE, row.names=1)

# 2. Correlations:
corrs = as.matrix(read.csv('./Results/2_Default_Parameters/metrics/average_correlations.csv', header=TRUE, row.names=1))
corrs_samples = read.csv('./Results/2_Default_Parameters/metrics/samples_correlations.csv', header=TRUE, row.names=1)
corrs_cellTypes = read.csv('./Results/2_Default_Parameters/metrics/cellTypes_correlations.csv', header=TRUE, row.names=1)





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
                        col = circlize::colorRamp2(c(0.1, 0.3, 0.5),
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
                              col = circlize::colorRamp2(c(0.1, 0.3, 0.5),
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
results_dir = './Results/2_Default_Parameters/proportions'
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
# 4.3.2. Read ground-truth proportions:
ground_truth = read.csv('./Data/bulk/CRC/NICs/data/new_metadata/cell_proportions.csv', row.names=1)
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




# -----------------------------------------------------
# --- COMPARE PREDICTION ABILITY TO NUMBER OF CELLS ---
# -----------------------------------------------------

# 0. Get cell counts info:
cell_counts = read.csv('./Data/bulk/CRC/NICs/data/new_metadata/cell_counts.csv', row.names=1)
total_cell_counts = colSums(cell_counts)

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

# 2. Add worst samples for each model:
df_rmse$label = df_rmse$samples
df_corr$label = df_corr$samples
for(model in unique(df_rmse$models)){
  mrmse = df_rmse[df_rmse$models==model, ]
  mcorr = df_corr[df_corr$models==model, ]
  
  rmse_samps = mrmse$samples[mrmse$values>0.2]
  corr_samps = mcorr$samples[mcorr$values<0.5]
  samps = intersect(rmse_samps, corr_samps)
  
  df_rmse$label[!df_rmse$samples%in%samps & df_rmse$models==model] = ''
  df_corr$label[!df_corr$samples%in%samps & df_corr$models==model] = ''
}


# 2. Plot how samples rmses correlate to counts, separated by method: (focus on the 5 best models)
models = c('CIBERSORTx', 'DigitalDLSorter', 'Scaden', 'AutoGeneS_nusvr', 'DWLS')#unique(df_rmse$models))
rmses_plots = list()
for(comparison_type in c('total_counts', 'cancer_counts', 'stromal_counts',
                         'immune_counts', 'others_counts')){
  fited_res = c()
  min_rmse = c()
  min_counts = c()
  for(model in models){
    fited_res = c(fited_res, paste('p =',
                                   round(cor(df_rmse$values[df_rmse$models==model],
                                       df_rmse[df_rmse$models==model, comparison_type]), 3)))
    min_rmse = c(min_rmse, min(df_rmse[df_rmse$models==model,'values']))
    min_counts = c(min_counts,
                   max(df_rmse[df_rmse$models==model, comparison_type]) -
                     (max(df_rmse[df_rmse$models==model, comparison_type])/10))
  }
  names(fited_res) = models
  
  rmses_plots[[comparison_type]] = ggplot2::ggplot(df_rmse[df_rmse$models%in%models,],
                                    ggplot2::aes_string(comparison_type, 'values', colour='models')) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(ggplot2::vars(models), scales='free') +
    ggrepel::geom_text_repel(ggplot2::aes(label=label), colour='black') +
    ggplot2::scale_colour_manual(values=methods_colours) +
    ggplot2::geom_hline(yintercept=0.2, linetype='dashed', colour='darkgrey') +
    ggplot2::geom_smooth(method='lm', colour='#80b3ff') +
    ggplot2::geom_text(data=data.frame(x=min_counts, y=min_rmse, lab=fited_res, models=names(fited_res)),
                       ggplot2::aes(x, y, label=lab), parse=FALSE, color='black', size=3) +
    ggplot2::theme(legend.position = 'none') + ggplot2::ylab('RMSE') + ggplot2::xlab('Cell counts')
}
#

# 3. Plot how samples correlations correlate to counts, separated by method:
models = c('CIBERSORTx', 'DigitalDLSorter', 'Scaden', 'AutoGeneS_nusvr', 'DWLS')#unique(df_rmse$models))
corrs_plots = list()
for(comparison_type in c('total_counts', 'cancer_counts', 'stromal_counts',
                         'immune_counts', 'others_counts')){
  fited_res = c()
  max_corr = c()
  min_counts = c()
  for(model in models){
    fited_res = c(fited_res, paste('p =',
                                   round(cor(df_corr$values[df_corr$models==model],
                                             df_corr[df_corr$models==model, comparison_type]), 3)))
    max_corr = c(max_corr, max(df_corr[df_corr$models==model,'values']))
    min_counts = c(min_counts,
                   max(df_corr[df_corr$models==model, comparison_type]) -
                     (max(df_corr[df_corr$models==model, comparison_type])/10))
  }
  names(fited_res) = models
  
  corrs_plots[[comparison_type]] = ggplot2::ggplot(df_corr[df_corr$models%in%models,],
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

# 0. Get cell counts info:
bulk_data = read.csv('./Data/bulk/CRC/NICs/data/data.csv', row.names=1)
total_read_counts = colSums(bulk_data)

# 1. Base dataframes for plots:
# 1.1. RMSEs
df_rmse_rc = rmses_samples[rmses_samples$cols=='uncorrected_gd',]
df_rmse_rc$read_counts = total_read_counts

# 1.2. Correlations
df_corr_rc = corrs_samples[corrs_samples$cols=='uncorrected_gd',]
df_corr_rc$read_counts = total_read_counts


# 2. Add worst samples for each model:
df_rmse_rc$label = df_rmse_rc$samples
df_corr_rc$label = df_corr_rc$samples
for(model in unique(df_rmse_rc$models)){
  mrmse = df_rmse_rc[df_rmse_rc$models==model, ]
  mcorr = df_corr_rc[df_corr_rc$models==model, ]
  
  rmse_samps = mrmse$samples[mrmse$values>0.2]
  corr_samps = mcorr$samples[mcorr$values<0.5]
  samps = intersect(rmse_samps, corr_samps)
  
  df_rmse_rc$label[!df_rmse_rc$samples%in%samps & df_rmse_rc$models==model] = ''
  df_corr_rc$label[!df_corr_rc$samples%in%samps & df_corr_rc$models==model] = ''
}


# 3. Plot how samples rmses correlate to counts, separated by method: (focus on the 5 best models)
models = c('CIBERSORTx', 'DigitalDLSorter', 'Scaden', 'AutoGeneS_nusvr', 'DWLS')#unique(df_rmse_rc$models))
fited_res = c()
min_rmse = c()
min_counts = c()
for(model in models){
  fited_res = c(fited_res, paste('p =',
                                 round(cor(df_rmse_rc$values[df_rmse_rc$models==model],
                                           df_rmse_rc[df_rmse_rc$models==model, 'read_counts']), 3)))
  min_rmse = c(min_rmse, min(df_rmse_rc[df_rmse_rc$models==model,'values']))
  min_counts = c(min_counts,
                 max(df_rmse_rc[df_rmse_rc$models==model, 'read_counts']) -
                   (max(df_rmse_rc[df_rmse_rc$models==model, 'read_counts'])/10))
}
names(fited_res) = models

ggplot2::ggplot(df_rmse_rc[df_rmse_rc$models%in%models,],
                ggplot2::aes_string('read_counts', 'values', colour='models')) +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(ggplot2::vars(models), scales='free') +
  ggrepel::geom_text_repel(ggplot2::aes(label=label), colour='black') +
  ggplot2::scale_colour_manual(values=methods_colours) +
  ggplot2::geom_hline(yintercept=0.2, linetype='dashed', colour='darkgrey') +
  ggplot2::geom_smooth(method='lm', colour='#80b3ff') +
  ggplot2::geom_text(data=data.frame(x=min_counts, y=min_rmse, lab=fited_res, models=names(fited_res)),
                     ggplot2::aes(x, y, label=lab), parse=FALSE, color='black', size=3) +
  ggplot2::theme(legend.position = 'none') + ggplot2::ylab('RMSE') + ggplot2::xlab('Cell counts')
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

# 0. Get cell counts info:
cell_props = read.csv('./Data/bulk/CRC/NICs/data/new_metadata/cell_proportions.csv', row.names=1)

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

# 2. Add worst samples for each model:
df_rmse_p$label = df_rmse_p$samples
df_corr_p$label = df_rmse_p$samples
for(model in unique(df_rmse_p$models)){
  mrmse = df_rmse_p[df_rmse_p$models==model, ]
  mcorr = df_corr_p[df_corr_p$models==model, ]
  
  rmse_samps = mrmse$samples[mrmse$values>0.2]
  corr_samps = mcorr$samples[mcorr$values<0.5]
  samps = intersect(rmse_samps, corr_samps)
  
  df_rmse_p$label[!df_rmse_p$samples%in%samps & df_rmse_p$models==model] = ''
  df_corr_p$label[!df_corr_p$samples%in%samps & df_corr_p$models==model] = ''
}


# 2. Plot how samples rmses correlate to counts, separated by method: (focus on the 5 best models)
models = c('CIBERSORTx', 'DigitalDLSorter', 'Scaden', 'AutoGeneS_nusvr', 'DWLS') #unique(df_rmse$models))
rmses_plots = list()
for(comparison_type in c('cancer_props', 'stromal_props', 'immune_props', 'others_props')){
  fited_res = c()
  min_rmse = c()
  min_counts = c()
  for(model in models){
    fited_res = c(fited_res, paste('p =',
                                   round(cor(df_rmse_p$values[df_rmse_p$models==model],
                                             df_rmse_p[df_rmse_p$models==model, comparison_type]), 3)))
    min_rmse = c(min_rmse, min(df_rmse_p[df_rmse_p$models==model,'values']))
    min_counts = c(min_counts,
                   max(df_rmse_p[df_rmse_p$models==model, comparison_type]) -
                     (max(df_rmse_p[df_rmse_p$models==model, comparison_type])/10))
  }
  names(fited_res) = models
  
  rmses_plots[[comparison_type]] = ggplot2::ggplot(df_rmse_p[df_rmse_p$models%in%models,],
                                                   ggplot2::aes_string(comparison_type, 'values', colour='models')) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(ggplot2::vars(models), scales='free') +
    ggrepel::geom_text_repel(ggplot2::aes(label=label), colour='black') +
    ggplot2::scale_colour_manual(values=methods_colours) +
    ggplot2::geom_hline(yintercept=0.2, linetype='dashed', colour='darkgrey') +
    ggplot2::geom_smooth(method='lm', colour='#80b3ff') +
    ggplot2::geom_text(data=data.frame(x=min_counts, y=min_rmse, lab=fited_res, models=names(fited_res)),
                       ggplot2::aes(x, y, label=lab), parse=FALSE, color='black', size=3) +
    ggplot2::theme(legend.position = 'none') + ggplot2::ylab('RMSE') + ggplot2::xlab('Cell counts')
}
#

# 3. Plot how samples correlations correlate to counts, separated by method:
models = c('CIBERSORTx', 'DigitalDLSorter', 'Scaden', 'AutoGeneS_nusvr', 'DWLS')#unique(df_corr_p$models))
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

