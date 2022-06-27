cell_type_colors = c('#011627', '#8b4513', '#f1558e', '#a5be00', '#800080', '#ff9f1c', '#2ec4b6', '#e71d36', '#3e673c', '#f9c22e', '#c4b1ae')
names(cell_type_colors) = c('Cancer cells', 'Stromal cells', 'Anti-Inflammatory macro/mono', 'Pro-Inflammatory macro/mono', 'Bcells',
                            'CD4 Tcells', 'Regulatory CD4 Tcells', 'CD8 Tcells', 'Proliferative Tcells', 'NK cells', 'Other cells')





# ----------------------------------------------
# --- HELPER PLOT FUNCTIONS OF TRAIN RESULTS ---
# ----------------------------------------------

# 1. Boxplots of the different values of a parameter for a method:
boxplots_parameter_method = function(df, parameter, fill_legend, ylim=c(0,2)){
  ggplot2::ggplot(df, ggplot2::aes_string(parameter, 'rmse', fill=parameter)) +
    ggplot2::geom_boxplot() + ggplot2::labs(x=fill_legend, y='RMSE')
    ggplot2::ylim(ylim[1], ylim[2]) + 
    ggplot2::theme(legend.position = 'none',
                   panel.background=ggplot2::element_rect(fill='white', colour='white'),
                   panel.grid.major=ggplot2::element_line(colour='grey'), panel.grid.minor=ggplot2::element_line(colour='grey')) +
    ggplot2::scale_fill_brewer(palette='Reds')
  
}

# 2. Scatter plot of rmses vs correlations, by cell-type, coloured by correction type
scatterPlot_rmsecorr_correction_celltype = function(proportions_list, ground_truth, train_samples,
                                                    legend_pos=c(.9,.15), ylim=c(-.11,1), xlim=c(0,.6)){
  cor_cts_methods = c()
  rmse_cts_methods = c()
  corrections_df = c()
  cellTypes_df = c()
  for(correction in names(proportions_list)){
    for(ct in cell_types){
      corrections_df = c(corrections_df, correction)
      cellTypes_df = c(cellTypes_df, ct)
      cor_cts_methods = c(cor_cts_methods, cor(as.numeric(ground_truth[ct,train_samples]),
                                                         as.numeric(proportions_list[[correction]][ct,train_samples]),
                                                         method = 'pearson'))
      rmse_cts_methods = c(rmse_cts_methods,Metrics::rmse(as.numeric(ground_truth[ct,train_samples]),
                                                          as.numeric(proportions_list[[correction]][ct,train_samples])))
    }
  }
  df_scatter_metrics_cellTypes = data.frame(Correction.Type=corrections_df, RMSE=rmse_cts_methods, Correlation=cor_cts_methods,
                                                 CellTypes=cellTypes_df)
  df_scatter_metrics_cellTypes[is.na(df_scatter_metrics_cellTypes)] = 0
  df_scatter_metrics_cellTypes$Correction.Type = factor(df_scatter_metrics_cellTypes$Correction.Type,
                                                        levels=c('uncorrected', 'corrected.Before', 'corrected.After'))
  ggplot2::ggplot(df_scatter_metrics_cellTypes, ggplot2::aes_string('RMSE', 'Correlation', colour='Correction.Type')) + 
    ggplot2::facet_wrap(ggplot2::vars(CellTypes)) +
    ggplot2::geom_point(size=3) + ggplot2::ylim(ylim) + ggplot2::xlim(xlim) +
    ggplot2::scale_colour_manual(values=c('#b30000', '#52854c', '#4e84c4')) +
    ggplot2::theme(legend.position = legend_pos)
}

# 3. Correction type plot:
scatterPlot_proportions_correction_type = function(proportions_list, ground_truth, cell_types = names(cell_type_colors),
                                                   lev_g=c('uncorrected', 'corrected.Before', 'corrected.After'),
                                                   cols_vector=cell_type_colors, plot_title='', xlim=c(0, 1),
                                                   ylim=c(0, 1), legend_pos=c(.87,.2), ncol=2, nrow=2, perCT=FALSE){
  
  # Data.frame for ggplot2:
  est_props = c()
  ground_vec = c()
  corrections = c()
  cell_type = c()
  for(correction in names(proportions_list)){
    for(samp in colnames(proportions_list[[correction]])){
      est_props = c(est_props, as.numeric(proportions_list[[correction]][cell_types, samp]))
      ground_vec = c(ground_vec, ground_truth[cell_types, samp])
      corrections = c(corrections, rep(correction, length(ground_truth[cell_types, samp])))
    }
    cell_type = c(cell_type, rep(cell_types, dim(proportions_list[[correction]])[2]))
  }
  dv = data.frame(cell_type = cell_type,
                  estimated_proportions = est_props,
                  ground_truth = ground_vec,
                  corrections = corrections)
  dv$corrections = factor(dv$corrections, levels=lev_g)
  dv$cell_type = factor(dv$cell_type, levels=cell_types)
  
  # First plot:
  plt = ggplot2::ggplot(data=dv, mapping=ggplot2::aes(ground_truth, estimated_proportions)) +
    ggplot2::geom_point(size=2, alpha=0.6, ggplot2::aes(colour=cell_type)) +
    ggplot2::geom_abline(intercept=0, slope=1, color="#404040", linetype='dashed') +
    ggplot2::xlim(xlim[1], xlim[2]) + ggplot2::ylim(ylim[1], ylim[2]) + ggplot2::ggtitle(plot_title)
  
  # Add facets, correlation and rmse according to perCT:
  if(perCT){
    plt = plt + ggplot2::facet_grid(ggplot2::vars(cell_type), ggplot2::vars(corrections))
    # plt = plt + ggpubr::stat_cor(ggplot2::aes(label= ..r.label..))
    corrs_anotate = c()
    rmses_anotate = c()
    corrections_anotate = c()
    cell_types_annotate = c()
    for(correction in names(proportions_list)){
      for(ct in cell_types){
        corr_val = cor(dv$ground_truth[dv$corrections==correction & dv$cell_type==ct],
                       dv$estimated_proportions[dv$corrections==correction & dv$cell_type==ct], method='pearson')
        rmse_val = Metrics::rmse(dv$ground_truth[dv$corrections==correction & dv$cell_type==ct],
                                 dv$estimated_proportions[dv$corrections==correction & dv$cell_type==ct])
        corr_anotate = paste("r: ", as.character(round(corr_val, 3)), sep='')
        rmse_anotate = paste('RMSE: ', as.character(round(rmse_val, 3)), sep='')
        corrs_anotate = c(corrs_anotate, corr_anotate)
        rmses_anotate = c(rmses_anotate, rmse_anotate)
        corrections_anotate = c(corrections_anotate, correction)
        cell_types_annotate = c(cell_types_annotate, ct)
      }
    }
    text_annotations = paste(corrs_anotate, rmses_anotate, sep=', ')
    dat_text = data.frame(label=text_annotations, corrections=corrections_anotate, cell_type=cell_types_annotate, stringsAsFactors=F)
    dat_text$corrections = factor(dat_text$corrections, levels=lev_g)
    dat_text$cell_type = factor(dat_text$cell_type, levels=cell_types)
    
    legend_pos = 'none'
  }
  else{
    plt = plt + ggplot2::facet_wrap(.~corrections, ncol=ncol, nrow=nrow)
    
    corrs_anotate = c()
    rmses_anotate = c()
    corrections_anotate = c()
    for(correction in names(proportions_list)){
      corr_val = cor(dv$ground_truth[dv$corrections==correction], dv$estimated_proportions[dv$corrections==correction], method='pearson')
      rmse_val = Metrics::rmse(dv$ground_truth[dv$corrections==correction], dv$estimated_proportions[dv$corrections==correction])
      corr_anotate = paste("r: ", as.character(round(corr_val,3 )), sep='')
      rmse_anotate = paste('RMSE: ', as.character(round(rmse_val,3 )), sep='')
      corrs_anotate = c(corrs_anotate, corr_anotate)
      rmses_anotate = c(rmses_anotate, rmse_anotate)
      corrections_anotate = c(corrections_anotate, correction)
    }
    text_annotations = paste(corrs_anotate, rmses_anotate, sep=', ')
    dat_text = data.frame(label=text_annotations, corrections=corrections_anotate, stringsAsFactors=F)
    dat_text$corrections = factor(dat_text$corrections, levels=lev_g)
    
    legend_pos = legend_pos
  }
  
  # Add correlation and rmse annotation:
  plt = plt + ggplot2::geom_text(data=dat_text, mapping=ggplot2::aes(x=xlim[2]/3.5, y=ylim[2]-0.01, label=label), size=3.2,
                                 parse=F, inherit.aes=FALSE)

  # Finalise plot:
  plt = plt + ggplot2::geom_smooth(method='lm', colour='#80b3ff') +
    ggplot2::scale_colour_manual(values=cols_vector) +
    ggplot2::labs(x='Ground Truth', y='Estimated Proportions', colour='Cell Types') +
    ggplot2::theme(panel.background=ggplot2::element_rect(fill='white', colour='white'),
                   panel.grid.major=ggplot2::element_line(colour='grey'), panel.grid.minor=ggplot2::element_line(colour='grey'),
                   legend.position = legend_pos, legend.background=ggplot2::element_rect(fill='white', size=.6,
                                                                                        linetype='solid', colour='black'),
                   legend.title=ggplot2::element_text(face='bold'), legend.text=ggplot2::element_text(size=9),
                   strip.text.y = ggplot2::element_blank(), strip.text.x = ggplot2::element_text(face='bold'))
  
  plt
}





# ----------------------------------------------
# --- HELPER PLOT FUNCTIONS OF TEST RESULTS ---
# ----------------------------------------------

# 1. Scatter plot of estimated vs ground-truth proportions
plot_estimated_vs_groundTruth_single_method = function(proportions, ground_truth, cell_types, cols_vector=cell_type_colors, annotate=TRUE,
                                                       plot_title='', xlim=c(0, 1), ylim=c(0, 1)){
  # Estimation + ground truth data.frame:
  est_props = c()
  ground_vec = c()
  for(samp in colnames(proportions)){
    est_props = c(est_props, as.numeric(proportions[cell_types, samp]))
    ground_vec = c(ground_vec, ground_truth[cell_types, samp])
  }
  dv = data.frame(cell_type = rep(cell_types, dim(proportions)[2]),
                  estimated_proportions = est_props,
                  ground_truth = ground_vec)

  # Plot:
  plt = ggplot2::ggplot(data=dv, mapping=ggplot2::aes(ground_truth, estimated_proportions)) +
    ggplot2::geom_point(size=2, alpha=0.6, ggplot2::aes(colour=cell_type)) +
    ggplot2::geom_abline(intercept=0, slope=1, color="#404040", linetype='dashed')
  
  # If annotation with correlation and RMSE metrics is wanted:
  if(annotate){
    corr_val = cor(dv$ground_truth, dv$estimated_proportions, method='pearson')
    rmse_val = Metrics::rmse(dv$ground_truth, dv$estimated_proportions)
    anotate = paste("r: ", as.character(round(corr_val,3 )), ', ', 'RMSE: ', as.character(round(rmse_val,3 )), sep='')
    
    plt = plt + ggplot2::annotate('text', x=xlim[2]/4, y=ylim[2]-0.01, label=anotate, parse=F)
  }
  
  plt = plt +
    ggplot2::geom_smooth(method='lm', colour='#80b3ff') +
    ggplot2::xlim(xlim[1], xlim[2]) + ggplot2::ylim(ylim[1], ylim[2]) + ggplot2::ggtitle(plot_title) +
    ggplot2::scale_colour_manual(values=cols_vector) +
    ggplot2::xlab('Ground Truth') + ggplot2::ylab('Estimated Proportions') + ggplot2::guides(color=ggplot2::guide_legend(title='Cell Types'))
  
  plt
}



# 2. Scatter plot of estimated vs ground-truth proportions by cell-type
plot_estimated_vs_groundTruth_perMethod_singleCT = function(proportions, ground_truth, cell_type, cols_vector,
                                                            plot_title=''){
  # Estimation + ground truth data.frame:
  est_props = c()
  ground_vec = c()
  methods_df = c()
  ymax = c()
  xmax = c()
  for(method in names(proportions)){
    samps = colnames(proportions[[method]])
    est_props = c(est_props, as.numeric(proportions[[method]][cell_type, samps]))
    ground_vec = c(ground_vec, as.numeric(ground_truth[cell_type, samps]))
    methods_df = c(methods_df, rep(method, length(samps)))
    ymax = c(ymax, max(as.numeric(proportions[[method]][cell_type, samps])))
    xrange = (max(as.numeric(ground_truth[cell_type, samps]))-
               min(as.numeric(ground_truth[cell_type, samps]))) / 2
    xmax = c(xmax, min(as.numeric(ground_truth[cell_type, samps])) + xrange)
  }
  dv = data.frame(methods = methods_df,
                  estimated_proportions = est_props,
                  ground_truth = ground_vec)
  # Metrics:
  corrs_anotate = c()
  rmses_anotate = c()
  cell_types_anotate = c()
  lm_fit_intercepts = c()
  lm_fit_slopes = c()
  for(method in names(proportions)){
    corr_val = cor(dv$ground_truth[dv$methods==method], dv$estimated_proportions[dv$methods==method], method='pearson')
    rmse_val = Metrics::rmse(dv$ground_truth[dv$methods==method], dv$estimated_proportions[dv$methods==method])
    corr_anotate = paste("p: ", as.character(round(corr_val,3 )), sep='')
    rmse_anotate = paste('RMSE: ', as.character(round(rmse_val,3 )), sep='')
    corrs_anotate = c(corrs_anotate, corr_anotate)
    rmses_anotate = c(rmses_anotate, rmse_anotate)
    cell_types_anotate = c(cell_types_anotate, method)
  }
  
  # Plot:
  plt = ggplot2::ggplot(data=dv, mapping=ggplot2::aes(ground_truth, estimated_proportions)) +
    ggplot2::geom_point(size=2, alpha=0.6, ggplot2::aes(colour=methods)) +
    ggplot2::geom_abline(intercept=0, slope=1, color="#404040", linetype='dashed') +
    ggplot2::geom_smooth(method='lm', colour='#80b3ff') +
    #ggplot2::xlim(xlim[1], xlim[2]) + ggplot2::ylim(ylim[1], ylim[2]) +
    ggplot2::ggtitle(plot_title) + ggplot2::facet_wrap(ggplot2::vars(methods), scales='free')
  
  
  
  corrs_dat_text = data.frame(label=paste(rmses_anotate, corrs_anotate, sep=' | '),
                              methods=cell_types_anotate, x=xmax, y=ymax, stringsAsFactors=F)
  plt = plt + ggplot2::geom_text(data=corrs_dat_text, mapping=ggplot2::aes(x=x, y=y, label=label),
                                 parse=F, inherit.aes=FALSE)
  #rmses_dat_text = data.frame(label=rmses_anotate, methods=cell_types_anotate, 
  #                            x=xmax, y=ymax-(ymax/1000),
  #                            stringsAsFactors=F)
  #plt = plt + ggplot2::geom_text(data=rmses_dat_text, mapping=ggplot2::aes(x=x, y=x, label=label),
  #                               parse=T, inherit.aes=FALSE)
  
  plt = plt + ggplot2::scale_colour_manual(values=cols_vector) + ggplot2::theme(legend.position='none') + 
    ggplot2::xlab('Ground Truth') + ggplot2::ylab('Estimated Proportions')
  
  plt
}