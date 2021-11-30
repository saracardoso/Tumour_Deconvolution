code_dir = './code/R'
for(file in list.files(code_dir, full.names=TRUE)) source(file)

# methods = c('AutoGeneS', 'BisqueRNA', 'BSeqSC', 'CIBERSORTx', 'DigitalDLSorter', 'DWLS', 'MOMF', 'MuSiC_woGrouping', 'MuSiC_wGrouping',
#             'Scaden', 'SCDC')
methods = c('AutoGeneS', 'BisqueRNA', 'BSeqSC', 'CIBERSORTx', 'DWLS', 'MOMF', 'MuSiC_woGrouping', 'MuSiC_wGrouping',
            'Scaden', 'SCDC')
cell_types = c('Cancer cells', 'Stromal cells', 'Anti-Inflammatory macro/mono', 'Pro-Inflammatory macro/mono', 'Bcells', 'CD4 Tcells',
               'Regulatory CD4 Tcells', 'CD8 Tcells', 'Proliferative Tcells', 'NK cells', 'Other cells')





# ---------------------
# --- TRAIN RESULTS ---
# ---------------------

train_results = list()
train_rmses_average_methods = c()
train_rmses_sd_methods = c()
train_rmses_samples_methods = c()
for(method in methods){
  if(method == 'DigitalDLSorter') next
  train_results[[method]] = read.csv(paste('./Results/1_Parameter_Optimization/train/', method, '.csv', sep=''), row.names=1)
  col_best = which.min(train_results[[method]]['MEAN',])
  nsamps = 16
  train_rmses_average_methods = c(train_rmses_average_methods, train_results[[method]]['MEAN',col_best])
  train_rmses_sd_methods = c(train_rmses_sd_methods, sd(train_results[[method]][1:nsamps, col_best]))
  train_rmses_samples_methods = rbind(train_rmses_samples_methods, train_results[[method]][1:nsamps, col_best])
}
names(train_rmses_average_methods) = methods
names(train_rmses_sd_methods) = methods
rownames(train_rmses_samples_methods) = methods
colnames(train_rmses_samples_methods) = rownames(train_results[[method]])[1:nsamps]





# --------------------
# --- TEST RESULTS ---
# --------------------

test_results_dir = './Results/1_Parameter_Optimization/test/proportions'
ground_truth_file = './Data/bulk/CRC/NICs/data/cell_proportions.csv'

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
test_rmses_sd_methods = c()
for(method in methods){
  test_rmses_full = rmse_dataset(predicted_proportions[[method]], ground_truth_proportions, cell_types)
  test_rmses_samples_methods = rbind(test_rmses_samples_methods, test_rmses_full$samples)
  test_rmses_average_methods = c(test_rmses_average_methods, test_rmses_full$dataset_average)
  test_rmses_sd_methods = c(test_rmses_sd_methods, sd(test_rmses_full$samples))
}
names(test_rmses_average_methods) = methods
names(test_rmses_sd_methods) = methods
rownames(test_rmses_samples_methods) = methods





# ---------------------
# --- TRAIN VS TEST ---
# ---------------------

# 1. Plot test vs train rmse averages, coloured by method:
df_averages = data.frame(rmse = c(train_rmses_average_methods, test_rmses_average_methods),
                         sd = c(train_rmses_sd_methods, test_rmses_sd_methods),
                         tt=c(rep('Train', length(methods)), rep('Test', length(methods))),
                         method = c(methods, methods))
df_averages$tt = factor(df_averages$tt, levels = c('Train', 'Test'))

ggplot2::ggplot(df_averages, ggplot2::aes(method, rmse, fill=tt)) +
  ggplot2::geom_bar(stat='identity', position=ggplot2::position_dodge()) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin=rmse-sd, ymax=rmse+sd), width=.2, position=ggplot2::position_dodge(.9)) +
  ggplot2::labs(x='Deconvolution Methods', y='RMSE', fill='') + ggplot2::ylim(0, 0.3) +
  ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, vjust=0.69)) +
  ggplot2::scale_fill_brewer(palette='Blues')





# -------------
# --- TRAIN ---
# -------------

# ---
# - 1. Boxplots of rmse values by correction type, for each method:
# ---
samples_train_rmses = c()
samples_train_correction_type = c()
samples_train_methods = c()
nsamps = 16
for(method in c('BisqueRNA', 'BSeqSC', 'CIBERSORTx', 'DWLS', 'MOMF', 'MuSiC_wGrouping', 'Scaden', 'SCDC')){
  combos = colnames(train_results[[method]])
  correction_type = rep('uncorrected', dim(train_results[[method]])[2])
  correction_type[grep('corrected[.]Before', combos)] = 'Corrected.Before'
  correction_type[grep('corrected[.]After', combos)] = 'Corrected.After'
  for(samp in rownames(train_results[[method]])[1:nsamps]){
    samples_train_rmses = c(samples_train_rmses, as.numeric(train_results[[method]][samp,]))
    samples_train_correction_type = c(samples_train_correction_type, correction_type)
    samples_train_methods = c(samples_train_methods, rep(method, length(combos)))
  }
}
df_correction_type = data.frame(rmse=samples_train_rmses, correction_type=samples_train_correction_type, method=samples_train_methods)
df_correction_type$correction_type = factor(df_correction_type$correction_type,
                                            levels = c('uncorrected', 'Corrected.Before', 'Corrected.After'))

ggplot2::ggplot(df_correction_type, ggplot2::aes(correction_type, rmse, fill=correction_type)) +
  ggplot2::geom_boxplot() + ggplot2::facet_wrap(ggplot2::vars(method)) +
  ggplot2::labs(x='', y='RMSE', fill='Type of Correction') +
  ggplot2::theme(axis.text.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank(),
                 panel.background=ggplot2::element_rect(fill='white', colour='white'),
                 panel.grid.major=ggplot2::element_line(colour='grey'), panel.grid.minor=ggplot2::element_line(colour='grey')) +
  ggplot2::scale_fill_manual(values=c('#b30000', '#52854c', '#4e84c4'))



# ---
# - 2. What parameter was most crucial to each method then? Check these boxplots, for each method, separated by the different parameters...
# ---
nsamps = 16

# 2.1. AutoGeneS
# 2.1.1. Get options values
n_options = dim(train_results$AutoGeneS)[2]
autogenes_model = rep('nusvr', n_options)
autogenes_gene.markers = rep('Yes', n_options)
autogenes_mode = rep('Fixed', n_options)
autogenes_n.iterations = rep('100', n_options)
autogenes_weights = rep('Correlation', n_options)
autogenes_model[grep('A2', colnames(train_results$AutoGeneS))] = 'nnls'
autogenes_model[grep('A3', colnames(train_results$AutoGeneS))] = 'linear'
autogenes_gene.markers[grep('B2', colnames(train_results$AutoGeneS))] = 'No'
autogenes_mode[grep('C2', colnames(train_results$AutoGeneS))] = 'Standard'
autogenes_n.iterations[grep('D2', colnames(train_results$AutoGeneS))] = '500'
autogenes_n.iterations[grep('D3', colnames(train_results$AutoGeneS))] = '1000'
autogenes_weights[grep('F2', colnames(train_results$AutoGeneS))] = 'Correlation+Distance'
autogenes_weights[grep('F3', colnames(train_results$AutoGeneS))] = 'Distance'
# 2.1.2. Get RMSEs of each sample:
autogenes_rmses = c()
for(samp in 1:nsamps){
  autogenes_rmses = c(autogenes_rmses, as.numeric(train_results$AutoGeneS[samp,]))
}
# 2.1.3. Construct data.frame
autogenes_df = data.frame(rmse=autogenes_rmses, model=rep(autogenes_model, nsamps), gene.markers=rep(autogenes_gene.markers, nsamps),
                          mode=rep(autogenes_mode, nsamps), n.iterations=rep(autogenes_n.iterations, nsamps),
                          weights=rep(autogenes_weights, nsamps))
autogenes_df$n.iterations = factor(autogenes_df$n.iterations, levels = c('100', '500', '1000'))
# 2.1.4. Construct boxplots:
autogenes_bxplt_model = boxplots_parameter_method(autogenes_df, 'model', 'Model', ylim=c(0,1))
autogenes_bxplt_gene.markers = boxplots_parameter_method(autogenes_df, 'gene.markers', 'Gene Markers', ylim=c(0,1))
autogenes_bxplt_mode = boxplots_parameter_method(autogenes_df, 'mode', 'Mode', ylim=c(0,1))
autogenes_bxplt_n.iterations = boxplots_parameter_method(autogenes_df, 'n.iterations', 'Nº Iterations', ylim=c(0,1))
autogenes_bxplt_weights = boxplots_parameter_method(autogenes_df, 'weights', 'Weights', ylim=c(0,1))
ggpubr::ggarrange(autogenes_bxplt_model, autogenes_bxplt_gene.markers, autogenes_bxplt_mode, autogenes_bxplt_n.iterations,
                  autogenes_bxplt_weights, ncol=3, nrow=2)

# 2.2. BisqueRNA
# 2.2.1. Get options values
n_options = dim(train_results$BisqueRNA)[2]
bisquerna_gene.markers = rep('Yes', n_options)
bisquerna_old.cpm = rep('After Subsetting Genes', n_options)
bisquerna_gene.markers[grep('A2', colnames(train_results$BisqueRNA))] = 'No'
bisquerna_old.cpm[grep('B2', colnames(train_results$BisqueRNA))] = 'Before Subsetting Genes'
# 2.2.2. Get RMSEs of each sample:
bisquerna_rmses = c()
for(samp in 1:nsamps){
  bisquerna_rmses = c(bisquerna_rmses, as.numeric(train_results$BisqueRNA[samp,]))
}
# 2.2.3. Construct data.frame
bisquerna_df = data.frame(rmse=bisquerna_rmses, gene.markers=rep(bisquerna_gene.markers, nsamps), old.cpm=rep(bisquerna_old.cpm, nsamps))
# 2.2.4. Construct boxplots:
bisquerna_bxplt_gene.markers = boxplots_parameter_method(bisquerna_df, 'gene.markers', 'Gene Markers', ylim=c(0,1))
bisquerna_bxplt_old.cpm = boxplots_parameter_method(bisquerna_df, 'old.cpm', 'When to Convert to CPM', ylim=c(0,1))
ggpubr::ggarrange(bisquerna_bxplt_gene.markers, bisquerna_bxplt_old.cpm, ncol=2, nrow=1)

# 2.3. BSeqSC
# 2.3.1. Get options values
n_options = dim(train_results$BSeqSC)[2]
bseqsc_ct.scale = rep('Yes', n_options)
bseqsc_permutations = rep('0', n_options)
bseqsc_ct.scale[grep('A2', colnames(train_results$BSeqSC))] = 'No'
bseqsc_permutations[grep('B2', colnames(train_results$BSeqSC))] = '50'
bseqsc_permutations[grep('B3', colnames(train_results$BSeqSC))] = '100'
bseqsc_permutations[grep('B4', colnames(train_results$BSeqSC))] = '500'
bseqsc_permutations[grep('B5', colnames(train_results$BSeqSC))] = '1000'
# 2.3.2. Get RMSEs of each sample:
bseqsc_rmses = c()
for(samp in 1:nsamps){
  bseqsc_rmses = c(bseqsc_rmses, as.numeric(train_results$BSeqSC[samp,]))
}
# 2.3.3. Construct data.frame
bseqsc_df = data.frame(rmse=bseqsc_rmses, ct.scale=rep(bseqsc_ct.scale, nsamps), permutations=rep(bseqsc_permutations, nsamps))
bseqsc_df$permutations = factor(bseqsc_df$permutations, levels = c('0', '50', '100', '500', '1000'))
# 2.3.4. Construct boxplots:
bseqsc_bxplt_gene.markers = boxplots_parameter_method(bseqsc_df, 'ct.scale', 'Scale raw counts', ylim=c(0,1))
bseqsc_bxplt_permutations = boxplots_parameter_method(bseqsc_df, 'permutations', 'Nº Permutations', ylim=c(0,1))
ggpubr::ggarrange(bseqsc_bxplt_gene.markers, bseqsc_bxplt_permutations, ncol=2, nrow=1)

# 2.4. CIBERSORTx
# 2.4.1. Get options values
n_options = dim(train_results$CIBERSORTx)[2]
cibersortx_batch.correction = rep('Yes', n_options)
cibersortx_permutations = rep('0', n_options)
cibersortx_min.expression = rep('0', n_options)
cibersortx_batch.correction[grep('A2', colnames(train_results$CIBERSORTx))] = 'No'
cibersortx_permutations[grep('B2', colnames(train_results$CIBERSORTx))] = '50'
cibersortx_permutations[grep('B3', colnames(train_results$CIBERSORTx))] = '100'
cibersortx_permutations[grep('B4', colnames(train_results$CIBERSORTx))] = '500'
cibersortx_permutations[grep('B5', colnames(train_results$CIBERSORTx))] = '1000'
cibersortx_min.expression[grep('C2', colnames(train_results$CIBERSORTx))] = '0.5'
# 2.4.2. Get RMSEs of each sample:
cibersortx_rmses = c()
for(samp in 1:nsamps){
  cibersortx_rmses = c(cibersortx_rmses, as.numeric(train_results$CIBERSORTx[samp,]))
}
# 2.4.3. Construct data.frame
cibersortx_df = data.frame(rmse=cibersortx_rmses, batch.correction=rep(cibersortx_batch.correction, nsamps),
                           permutations=rep(cibersortx_permutations, nsamps), min.expression=rep(cibersortx_min.expression, nsamps))
cibersortx_df$permutations = factor(cibersortx_df$permutations, levels = c('0', '50', '100', '500', '1000'))
cibersortx_df$min.expression = factor(cibersortx_df$min.expression, levels = c('0', '0.5'))
# 2.4.4. Construct boxplots:
cibersortx_bxplt_batch.correction = boxplots_parameter_method(cibersortx_df, 'batch.correction', 'Batch Correction', ylim=c(0,1))
cibersortx_bxplt_permutations = boxplots_parameter_method(cibersortx_df, 'permutations', 'Nº Permutations', ylim=c(0,1))
cibersortx_bxplt_min.expression = boxplots_parameter_method(cibersortx_df, 'min.expression', 'Minimum Expression', ylim=c(0,1))
ggpubr::ggarrange(cibersortx_bxplt_batch.correction, cibersortx_bxplt_permutations, cibersortx_bxplt_min.expression, ncol=2, nrow=2)

# 2.5. DigitalDLSorter

# 2.6. DWLS
# 2.6.1. Get options values
n_options = dim(train_results$DWLS)[2]
dwls_diff.cutoff = rep('0.5', n_options)
dwls_pval.cutoff = rep('0.01', n_options)
dwls_diff.cutoff[grep('A2', colnames(train_results$DWLS))] = '0.8'
dwls_diff.cutoff[grep('A3', colnames(train_results$DWLS))] = '1'
dwls_pval.cutoff[grep('B2', colnames(train_results$DWLS))] = '0.05'
# 2.6.2. Get RMSEs of each sample:
dwls_rmses = c()
for(samp in 1:nsamps){
  dwls_rmses = c(dwls_rmses, as.numeric(train_results$DWLS[samp,]))
}
# 2.6.3. Construct data.frame
dwls_df = data.frame(rmse=dwls_rmses, diff.cutoff=rep(dwls_diff.cutoff, nsamps), pval.cutoff=rep(dwls_pval.cutoff, nsamps))
dwls_df$diff.cutoff = factor(dwls_df$diff.cutoff, levels = c('0.5', '0.8', '1'))
dwls_df$pval.cutoff = factor(dwls_df$pval.cutoff, levels = c('0.01', '0.05'))
# 2.6.4. Construct boxplots:
dwls_bxplt_diff.cutoff = boxplots_parameter_method(dwls_df, 'diff.cutoff', 'Minimum Fold Change', ylim=c(0,1))
dwls_bxplt_pval.cutoff = boxplots_parameter_method(dwls_df, 'pval.cutoff', 'p-Value', ylim=c(0,1))
ggpubr::ggarrange(dwls_bxplt_diff.cutoff, dwls_bxplt_pval.cutoff, ncol=2, nrow=1)

# 2.7. MOMF
# 2.7.1. Get options values
n_options = dim(train_results$MOMF)[2]
momf_method = rep('KL', n_options)
momf_rho = rep('1', n_options)
momf_num.iter = rep('100', n_options)
momf_method[grep('A2', colnames(train_results$MOMF))] = 'IS'
momf_rho[grep('B2', colnames(train_results$MOMF))] = '2'
momf_rho[grep('B3', colnames(train_results$MOMF))] = '5'
momf_num.iter[grep('C2', colnames(train_results$MOMF))] = '500'
momf_num.iter[grep('C3', colnames(train_results$MOMF))] = '1000'
# 2.7.2. Get RMSEs of each sample:
momf_rmses = c()
for(samp in 1:nsamps){
  momf_rmses = c(momf_rmses, as.numeric(train_results$MOMF[samp,]))
}
# 2.7.3. Construct data.frame
momf_df = data.frame(rmse=momf_rmses, method=rep(momf_method, nsamps), rho=rep(momf_rho, nsamps), num.iter=rep(momf_num.iter, nsamps))
momf_df$rho = factor(momf_df$rho, levels = c('1', '2', '5'))
momf_df$num.iter = factor(momf_df$num.iter, levels = c('100', '500', '1000'))
# 2.7.4. Construct boxplots:
momf_bxplt_method = boxplots_parameter_method(momf_df, 'method', 'Method', ylim=c(0,1))
momf_bxplt_rho = boxplots_parameter_method(momf_df, 'rho', 'rho', ylim=c(0,1))
momf_bxplt_num.iter = boxplots_parameter_method(momf_df, 'num.iter', 'Iterations', ylim=c(0,1))
ggpubr::ggarrange(momf_bxplt_method, momf_bxplt_rho, momf_bxplt_num.iter, ncol=2, nrow=2)

# 2.8. MuSiC_woGrouping
# 2.8.1. Get options values
n_options = dim(train_results$MuSiC_woGrouping)[2]
musicwo_n.iterations = rep('100', n_options)
musicwo_center = rep('Yes', n_options)
musicwo_normalise = rep('Yes', n_options)
musicwo_wo.method = rep('weighted', n_options)
musicwo_wo.gene.markers = rep('Yes', n_options)
musicwo_n.iterations[grep('A2', colnames(train_results$MuSiC_woGrouping))] = '500'
musicwo_n.iterations[grep('A3', colnames(train_results$MuSiC_woGrouping))] = '1000'
musicwo_center[grep('B2', colnames(train_results$MuSiC_woGrouping))] = 'No'
musicwo_normalise[grep('C2', colnames(train_results$MuSiC_woGrouping))] = 'No'
musicwo_wo.method[grep('D2', colnames(train_results$MuSiC_woGrouping))] = 'All Genes'
musicwo_wo.gene.markers[grep('E2', colnames(train_results$MuSiC_woGrouping))] = 'No'
# 2.8.2. Get RMSEs of each sample:
musicwo_rmses = c()
for(samp in 1:nsamps){
  musicwo_rmses = c(musicwo_rmses, as.numeric(train_results$MuSiC_woGrouping[samp,]))
}
# 2.8.3. Construct data.frame
musicwo_df = data.frame(rmse=musicwo_rmses, n.iterations=rep(musicwo_n.iterations, nsamps), center=rep(musicwo_center, nsamps),
                     normalise=rep(musicwo_normalise, nsamps), method=rep(musicwo_wo.method, nsamps),
                     gene.markers=rep(musicwo_wo.gene.markers, nsamps))
musicwo_df$n.iterations = factor(musicwo_df$n.iterations, levels = c('100', '500', '1000'))
musicwo_df$center = factor(musicwo_df$center, levels = c('Yes', 'No'))
musicwo_df$normalise = factor(musicwo_df$normalise, levels = c('Yes', 'No'))
musicwo_df$gene.markers = factor(musicwo_df$gene.markers, levels = c('Yes', 'No'))
# 2.8.4. Construct boxplots:
musicwo_bxplt_n.iterations = boxplots_parameter_method(musicwo_df, 'n.iterations', 'Iterations', ylim=c(0,1))
musicwo_bxplt_center = boxplots_parameter_method(musicwo_df, 'center', 'Center', ylim=c(0,1))
musicwo_bxplt_normalise = boxplots_parameter_method(musicwo_df, 'normalise', 'Normalise', ylim=c(0,1))
musicwo_bxplt_method = boxplots_parameter_method(musicwo_df, 'method', 'Method', ylim=c(0,1))
musicwo_bxplt_gene.markers = boxplots_parameter_method(musicwo_df, 'gene.markers', 'Gene Markers', ylim=c(0,1))
ggpubr::ggarrange(musicwo_bxplt_n.iterations, musicwo_bxplt_center, musicwo_bxplt_normalise, musicwo_bxplt_method,
                  musicwo_bxplt_gene.markers, ncol=3, nrow=2)

# 2.9. MuSiC_wGrouping
# 2.9.1. Get options values
n_options = dim(train_results$MuSiC_wGrouping)[2]
musicw_n.iterations = rep('100', n_options)
musicw_center = rep('Yes', n_options)
musicw_normalise = rep('Yes', n_options)
musicw_n.iterations[grep('A2', colnames(train_results$MuSiC_wGrouping))] = '500'
musicw_n.iterations[grep('A3', colnames(train_results$MuSiC_wGrouping))] = '1000'
musicw_center[grep('B2', colnames(train_results$MuSiC_wGrouping))] = 'No'
musicw_normalise[grep('C2', colnames(train_results$MuSiC_wGrouping))] = 'No'
# 2.9.2. Get RMSEs of each sample:
musicw_rmses = c()
for(samp in 1:nsamps){
  musicw_rmses = c(musicw_rmses, as.numeric(train_results$MuSiC_wGrouping[samp,]))
}
# 2.9.3. Construct data.frame
musicw_df = data.frame(rmse=musicw_rmses, n.iterations=rep(musicw_n.iterations, nsamps), center=rep(musicw_center, nsamps),
                        normalise=rep(musicw_normalise, nsamps))
musicw_df$n.iterations = factor(musicw_df$n.iterations, levels = c('100', '500', '1000'))
musicw_df$center = factor(musicw_df$center, levels = c('Yes', 'No'))
musicw_df$normalise = factor(musicw_df$normalise, levels = c('Yes', 'No'))
# 2.9.4. Construct boxplots:
musicw_bxplt_n.iterations = boxplots_parameter_method(musicw_df, 'n.iterations', 'Iterations', ylim=c(0,1))
musicw_bxplt_center = boxplots_parameter_method(musicw_df, 'center', 'Center', ylim=c(0,1))
musicw_bxplt_normalise = boxplots_parameter_method(musicw_df, 'normalise', 'Normalise', ylim=c(0,1))
ggpubr::ggarrange(musicw_bxplt_n.iterations, musicw_bxplt_center, musicw_bxplt_normalise, ncol=2, nrow=2)

# 2.10. Scaden
# 2.10.1. Get options values
n_options = dim(train_results$Scaden)[2]
scaden_min.expression = rep('0.1', n_options)
scaden_batch.size = rep('50', n_options)
scaden_learning.rate = rep('0.0001', n_options)
scaden_n.steps = rep('100', n_options)
scaden_min.expression[grep('A2', colnames(train_results$Scaden))] = '1'
scaden_min.expression[grep('A3', colnames(train_results$Scaden))] = '5'
scaden_batch.size[grep('B2', colnames(train_results$Scaden))] = '128'
scaden_batch.size[grep('B3', colnames(train_results$Scaden))] = '200'
scaden_learning.rate[grep('C2', colnames(train_results$Scaden))] = '0.001'
scaden_learning.rate[grep('C3', colnames(train_results$Scaden))] = '0.01'
scaden_n.steps[grep('D2', colnames(train_results$Scaden))] = '500'
scaden_n.steps[grep('D3', colnames(train_results$Scaden))] = '1000'
# 2.10.2. Get RMSEs of each sample:
scaden_rmses = c()
for(samp in 1:nsamps){
  scaden_rmses = c(scaden_rmses, as.numeric(train_results$Scaden[samp,]))
}
# 2.10.3. Construct data.frame
scaden_df = data.frame(rmse=scaden_rmses, min.expression=rep(scaden_min.expression, nsamps), batch.size=rep(scaden_batch.size, nsamps),
                       learning.rate=rep(scaden_learning.rate, nsamps), n.steps=rep(scaden_n.steps, nsamps))
scaden_df$min.expression = factor(scaden_df$min.expression, levels = c('0.1', '1', '5'))
scaden_df$batch.size = factor(scaden_df$batch.size, levels = c('50', '128', '200'))
scaden_df$learning.rate = factor(scaden_df$learning.rate, levels = c('0.0001', '0.001', '0.01'))
scaden_df$n.steps = factor(scaden_df$n.steps, levels = c('100', '500', '1000'))
# 2.10.4. Construct boxplots:
scaden_bxplt_min.expression = boxplots_parameter_method(scaden_df, 'min.expression', 'Minimum Expression', ylim=c(0,1))
scaden_bxplt_batch.size = boxplots_parameter_method(scaden_df, 'batch.size', 'Batch Size', ylim=c(0,1))
scaden_bxplt_learning.rate = boxplots_parameter_method(scaden_df, 'learning.rate', 'Learning Rate', ylim=c(0,1))
scaden_bxplt_n.steps = boxplots_parameter_method(scaden_df, 'n.steps', 'Nº Steps', ylim=c(0,1))
ggpubr::ggarrange(scaden_bxplt_min.expression, scaden_bxplt_batch.size, scaden_bxplt_learning.rate, scaden_bxplt_n.steps, ncol=2, nrow=2)

# 2.11. SCDC
# 2.11.1. Get options values
n_options = dim(train_results$Scaden)[2]
scdc_search.length = rep('0.01', n_options)
scdc_n.iterations = rep('100', n_options)
scdc_search.length[grep('A2', colnames(train_results$Scaden))] = '0.05'
scdc_search.length[grep('A3', colnames(train_results$Scaden))] = '0.1'
scdc_n.iterations[grep('B2', colnames(train_results$Scaden))] = '500'
scdc_n.iterations[grep('B3', colnames(train_results$Scaden))] = '1000'
# 2.11.2. Get RMSEs of each sample:
scdc_rmses = c()
for(samp in 1:nsamps){
  scdc_rmses = c(scdc_rmses, as.numeric(train_results$Scaden[samp,]))
}
# 2.11.3. Construct data.frame
scdc_df = data.frame(rmse=scdc_rmses, search.length=rep(scdc_search.length, nsamps), n.iterations=rep(scdc_n.iterations, nsamps))
scdc_df$search.length = factor(scdc_df$search.length, levels = c('0.01', '0.05', '0.1'))
scdc_df$n.iterations = factor(scdc_df$n.iterations, levels = c('100', '500', '1000'))
# 2.11.4. Construct boxplots:
scdc_bxplt_search.length = boxplots_parameter_method(scdc_df, 'search.length', 'Search Length', ylim=c(0,1))
scdc_bxplt_n.iterations = boxplots_parameter_method(scdc_df, 'n.iterations', 'Iterations', ylim=c(0,1))
ggpubr::ggarrange(scdc_bxplt_search.length, scdc_bxplt_n.iterations, ncol=2, nrow=1)



# ---
# - 3. Check proportions between correction types to see if the fact that we are not seeing better estimations overall with the corrections
# ---  happens across all cell-types or some cell-types are being indeed better estimated but it comes at the big expense of others (if so, it
# ---  could be a good reason to try out ensemble of methods)
# Train data:
train = readRDS('./Data/bulk/CRC/NICs/train_test_sets/train.Rdata')
train_bulk_file_scaden = './Data/bulk/CRC/NICs/train_test_sets/train.txt'
# References:
references = list()
references$uncorrected = read_scRNAseq_reference('./Data/CRC_reference/UB_matrix_CSV.csv', './Data/CRC_reference/metadata.csv')
invisible(gc())
references$corrected = read_scRNAseq_reference('./Data/CRC_reference/B_matrix_CSV.csv', './Data/CRC_reference/metadata.csv')
invisible(gc())
# Gene Markers:
markers =  jsonlite::read_json('./Data/CRC_reference/markers.json', simplifyVector=T)
# RNA bias correction vector
rna_correction_df = read.csv('./Data/CRC_reference/bias_values.csv', row.names=1)
rna_correction_vec = rna_correction_df[,1]
names(rna_correction_vec) = rownames(rna_correction_df)
# Names of the train samples
train_samples = c('S52', 'NIC4', 'NIC23', 'NIC24', 'NIC25', 'NIC12nova', 'NIC29', 'NIC15', 'NIC18', 'NIC17', 'NIC3', 'NIC22', 'NIC19',
                  'NIC20', 'NIC13', 'NIC7')


# 3.1. BisqueRNA
correction_type_results_BisqueRNA = list()
# Best combo: corrected.Before_A1B2
which.min(train_results$BisqueRNA['MEAN',])
# Best combo for no correction: uncorrected_A1B2
which.min(train_results$BisqueRNA['MEAN',grep('uncorrected', colnames(train_results$BisqueRNA))])
# Best combo for correction after estimation: corrected.After_A1B1
which.min(train_results$BisqueRNA['MEAN',grep('corrected.After', colnames(train_results$BisqueRNA))])
# 3.1.1. Run BisqueRNA for uncorrected_A1B2, corrected.Before_A1B2 and corrected.After_A1B1
# - uncorrected:
correction_type_results_BisqueRNA$uncorrected = Estimate.Proportions(train, references$uncorrected, 'Deconv_cellTypes', 'BisqueRNA',
                                                                     gene.markers=markers, use.overlap=FALSE, old.cpm=FALSE)
invisible(gc())
# - corrected before:
correction_type_results_BisqueRNA$corrected.Before = Estimate.Proportions(train, references$corrected, 'Deconv_cellTypes', 'BisqueRNA',
                                                                     gene.markers=markers, use.overlap=FALSE, old.cpm=FALSE)
invisible(gc())
# - corrected After:
res_temp = Estimate.Proportions(train, references$uncorrected, 'Deconv_cellTypes', 'BisqueRNA',
                                                                     gene.markers=markers, use.overlap=FALSE, old.cpm=TRUE)
correction_type_results_BisqueRNA$corrected.After = correct_fractions_mRNABias(res_temp, rna_correction_vec)
invisible(gc())
# 3.1.2. Save results:
saveRDS(correction_type_results_BisqueRNA, './Results/1_Parameter_Optimization/train/correction_type_analysis_res/BisqueRNA.Rdata')
# 3.1.3. Plot RMSE vs correlation by each cell-types, coloured by correction type:
scatterPlot_rmsecorr_correction_celltype(correction_type_results_BisqueRNA, ground_truth_proportions, train_samples,
                                         ylim=c(-.2,1), xlim=c(0,.65))


# 3.2. BSeqSC
correction_type_results_BSeqSC = list()
# Best combo: uncorrected_A2B1
which.min(train_results$BSeqSC['MEAN',])
# Best combo for correction before estimation: corrected.Before_A2B1
which.min(train_results$BSeqSC['MEAN',grep('corrected.Before', colnames(train_results$BSeqSC))])
# Best combo for correction after estimation: corrected.After_A2B1
which.min(train_results$BSeqSC['MEAN',grep('corrected.After', colnames(train_results$BSeqSC))])
# 3.2.1. Run BSeqSC for uncorrected_A2B1, corrected.Before_A2B1 and corrected.After_A2B1
# - uncorrected:
res_temp = t(read.csv('./Results/1_Parameter_Optimization/train/BSeqSC_web/uncorrected_A2B1.csv', row.names=1))
res_temp = res_temp[1:11,]
rownames(res_temp) = gsub('[_]', ' ', rownames(res_temp))
rownames(res_temp)[3:4] = c('Anti-Inflammatory macro/mono', 'Pro-Inflammatory macro/mono')
correction_type_results_BSeqSC$uncorrected = res_temp[cell_types,]
invisible(gc())
# - corrected before:
res_temp = t(read.csv('./Results/1_Parameter_Optimization/train/BSeqSC_web/corrected.Before_A2B1.csv', row.names=1))
res_temp = res_temp[1:11,]
rownames(res_temp) = gsub('[_]', ' ', rownames(res_temp))
rownames(res_temp)[3:4] = c('Anti-Inflammatory macro/mono', 'Pro-Inflammatory macro/mono')
correction_type_results_BSeqSC$corrected.Before = res_temp[cell_types,]
invisible(gc())
# - corrected After:
res_temp = t(read.csv('./Results/1_Parameter_Optimization/train/BSeqSC_web/uncorrected_A2B1.csv', row.names=1))
res_temp = res_temp[1:11,]
rownames(res_temp) = gsub('[_]', ' ', rownames(res_temp))
rownames(res_temp)[3:4] = c('Anti-Inflammatory macro/mono', 'Pro-Inflammatory macro/mono')
correction_type_results_BSeqSC$corrected.After = correct_fractions_mRNABias(res_temp, rna_correction_vec)
invisible(gc())
# 3.2.2. Save results:
saveRDS(correction_type_results_BSeqSC, './Results/1_Parameter_Optimization/train/correction_type_analysis_res/BSeqSC.Rdata')
# 3.2.3. Plot RMSE vs correlation by each cell-types, coloured by correction type:
scatterPlot_rmsecorr_correction_celltype(correction_type_results_BSeqSC, ground_truth_proportions, train_samples,
                                         ylim=c(-.5,1), xlim=c(0,.7))


# 3.3. CIBERSORTx
correction_type_results_CIBERSORTx = list()
# Best combo: uncorrected_A2B1C2
which.min(train_results$CIBERSORTx['MEAN',])
# Best combo for correction before estimation: corrected.Before_A2B1C2
which.min(train_results$CIBERSORTx['MEAN',grep('corrected.Before', colnames(train_results$CIBERSORTx))])
# Best combo for correction after estimation: corrected.After_A2B1C2
which.min(train_results$CIBERSORTx['MEAN',grep('corrected.After', colnames(train_results$CIBERSORTx))])
# 3.3.1. Run CIBERSORTx for uncorrected_A2B1C2, corrected.Before_A2B1C2 and corrected.After_A2B1C2
# - uncorrected:
res_temp = t(read.csv('./Results/1_Parameter_Optimization/train/CIBERSORTx_web/uncorrected_A2B1C2.csv', row.names=1))
res_temp = res_temp[1:11,]
rownames(res_temp) = gsub('[_]', ' ', rownames(res_temp))
rownames(res_temp)[3:4] = c('Anti-Inflammatory macro/mono', 'Pro-Inflammatory macro/mono')
correction_type_results_CIBERSORTx$uncorrected = res_temp[cell_types,]
invisible(gc())
# - corrected before:
res_temp = t(read.csv('./Results/1_Parameter_Optimization/train/CIBERSORTx_web/corrected.Before_A2B1C2.csv', row.names=1))
res_temp = res_temp[1:11,]
rownames(res_temp) = gsub('[_]', ' ', rownames(res_temp))
rownames(res_temp)[3:4] = c('Anti-Inflammatory macro/mono', 'Pro-Inflammatory macro/mono')
correction_type_results_CIBERSORTx$corrected.Before = res_temp[cell_types,]
invisible(gc())
# - corrected After:
res_temp = t(read.csv('./Results/1_Parameter_Optimization/train/CIBERSORTx_web/uncorrected_A2B1C2.csv', row.names=1))
res_temp = res_temp[1:11,]
rownames(res_temp) = gsub('[_]', ' ', rownames(res_temp))
rownames(res_temp)[3:4] = c('Anti-Inflammatory macro/mono', 'Pro-Inflammatory macro/mono')
correction_type_results_CIBERSORTx$corrected.After = correct_fractions_mRNABias(res_temp, rna_correction_vec)
invisible(gc())
# 3.3.2. Save results:
saveRDS(correction_type_results_CIBERSORTx, './Results/1_Parameter_Optimization/train/correction_type_analysis_res/CIBERSORTx.Rdata')
# 3.3.3. Plot RMSE vs correlation by each cell-types, coloured by correction type:
scatterPlot_rmsecorr_correction_celltype(correction_type_results_CIBERSORTx, ground_truth_proportions, train_samples,
                                         ylim=c(-.2,1), xlim=c(0,.4))


# 3.4. DWLS
correction_type_results_DWLS = list()
# Best combo: uncorrected_A3B1
which.min(train_results$DWLS['MEAN',])
# Best combo for correction after estimation: corrected.After_A3B2
which.min(train_results$DWLS['MEAN',grep('corrected.After', colnames(train_results$DWLS))])
# 3.4.1. Run DWLS for uncorrected_A3B1 and corrected.After_A3B2
# - uncorrected:
correction_type_results_DWLS$uncorrected = Estimate.Proportions(train, references$uncorrected, 'Deconv_cellTypes', 'DWLS',
                                                                1, 0.01, rna.bias=FALSE, bias.vec=NULL)
correction_type_results_DWLS$uncorrected[correction_type_results_DWLS$uncorrected<0] = 0
invisible(gc())
# - corrected After:
res_temp = Estimate.Proportions(train, references$uncorrected, 'Deconv_cellTypes', 'DWLS',
                                1, 0.05, rna.bias=FALSE, bias.vec=NULL)
correction_type_results_DWLS$corrected.After = correct_fractions_mRNABias(res_temp, rna_correction_vec)
invisible(gc())
# 3.4.2. Save results:
saveRDS(correction_type_results_DWLS, './Results/1_Parameter_Optimization/train/correction_type_analysis_res/DWLS.Rdata')
# 3.4.3. Plot RMSE vs correlation by each cell-types, coloured by correction type:
scatterPlot_rmsecorr_correction_celltype(correction_type_results_DWLS, ground_truth_proportions, train_samples,
                                         ylim=c(-.2,1), xlim=c(0,.4))


# 3.5. MOMF
correction_type_results_MOMF = list()
# Best combo: corrected.Before_A2B3C3
which.min(train_results$MOMF['MEAN',])
# Best combo for no correction: 
which.min(train_results$MOMF['MEAN',grep('uncorrected', colnames(train_results$MOMF))])
# Best combo for correction after estimation: 
which.min(train_results$MOMF['MEAN',grep('corrected.After', colnames(train_results$MOMF))])
# 3.5.1. Run MOMF for uncorrected_A1B3C2, corrected.Before_A2B3C3 and corrected.After_A1B3C2
# - uncorrected:
res_temp = Estimate.Proportions(train, references$uncorrected, 'Deconv_cellTypes', 'MOMF', 'KL',
                                                                5, 500, FALSE, NULL)
correction_type_results_MOMF$uncorrected = res_temp
invisible(gc())
# - corrected Before:
correction_type_results_MOMF$corrected.Before = Estimate.Proportions(train, references$corrected, 'Deconv_cellTypes', 'MOMF', 'IS',
                                                                5, 1000, FALSE, NULL)
invisible(gc())
# - corrected After:
correction_type_results_MOMF$corrected.After = correct_fractions_mRNABias(res_temp, rna_correction_vec)
invisible(gc())
# 3.5.2. Save results:
saveRDS(correction_type_results_MOMF, './Results/1_Parameter_Optimization/train/correction_type_analysis_res/MOMF.Rdata')
# 3.5.3. Plot RMSE vs correlation by each cell-types, coloured by correction type:
scatterPlot_rmsecorr_correction_celltype(correction_type_results_MOMF, ground_truth_proportions, train_samples,
                                         ylim=c(-.5,1), xlim=c(0,.65))


# 3.6. MuSiC_wGrouping
correction_type_results_MuSiC_wGrouping = list()
# Best combo: uncorrected_A1B2C1
which.min(train_results$MuSiC_wGrouping['MEAN',])
# Best combo for correction before estimation: corrected.Before_A1B2C2 
which.min(train_results$MuSiC_wGrouping['MEAN',grep('corrected.Before', colnames(train_results$MuSiC_wGrouping))])
# Best combo for correction after estimation: 
which.min(train_results$MuSiC_wGrouping['MEAN',grep('corrected.After', colnames(train_results$MuSiC_wGrouping))])
# 3.6.1. Run MuSiC_wGrouping for uncorrected_A1B2C1, corrected.Before_A1B2C2 and corrected.After_A1B2C1
clusters = list(C1=c("Cancer cells", "Stromal cells"), C2=c("Anti-Inflammatory macro/mono", "Pro-Inflammatory macro/mono", "Other cells"),
                C3=c("Bcells", "Regulatory CD4 Tcells", "CD4 Tcells", "Proliferative Tcells", "NK cells", "CD8 Tcells"))
# - uncorrected:
res_temp = Estimate.Proportions(train, references$uncorrected, 'Deconv_cellTypes', 'MuSiC', music.method='with.Grouping',
                                n.iterations=100, center=FALSE, normalize=TRUE,
                                w.celltypes.Cluster.List=clusters, w.markers.log2FC=.8,
                                w.markers.pval=.01, w.markers.only.Pos=TRUE)
correction_type_results_MuSiC_wGrouping$uncorrected = res_temp
invisible(gc())
# - corrected Before:
correction_type_results_MuSiC_wGrouping$corrected.Before = Estimate.Proportions(train, references$corrected, 'Deconv_cellTypes', 'MuSiC',
                                                                                music.method='with.Grouping',
                                                                                n.iterations=100, center=FALSE, normalize=FALSE,
                                                                                w.celltypes.Cluster.List=clusters, w.markers.log2FC=.8,
                                                                                w.markers.pval=.01, w.markers.only.Pos=TRUE)
invisible(gc())
# - corrected After:
correction_type_results_MuSiC_wGrouping$corrected.After = correct_fractions_mRNABias(res_temp, rna_correction_vec)
invisible(gc())
# 3.6.2. Save results:
saveRDS(correction_type_results_MuSiC_wGrouping,
        './Results/1_Parameter_Optimization/train/correction_type_analysis_res/MuSiC_wGrouping.Rdata')
# 3.6.3. Plot RMSE vs correlation by each cell-types, coloured by correction type:
scatterPlot_rmsecorr_correction_celltype(correction_type_results_MuSiC_wGrouping, ground_truth_proportions, train_samples,
                                         ylim=c(-.6,1), xlim=c(0,.7))


# 3.7. Scaden
correction_type_results_Scaden = list()
# Best combo: corrected.After_A1B1C2D3
which.min(train_results$Scaden['MEAN',])
# Best combo for uncorrectedn: uncorrected_A3B1C2D1 
which.min(train_results$Scaden['MEAN',grep('uncorrected', colnames(train_results$Scaden))])
# Best combo for correction before estimation: corrected.Before_A1B2C2D2
which.min(train_results$Scaden['MEAN',grep('corrected.Before', colnames(train_results$Scaden))])
# 3.7.1. Run Scaden for uncorrected_A1B2C1, corrected.Before_A1B2C2D2 and corrected.After_A1B1C2D3
cell_types_scaden = cell_types
names(cell_types_scaden) = cell_types
cell_types_scaden[3:4] = c('Anti Inflammatory macro mono', 'Pro Inflammatory macro mono')
rna_correction_vec_scaden = rna_correction_vec
names(rna_correction_vec_scaden) = cell_types_scaden[match(names(rna_correction_vec), names(cell_types_scaden))]
# - uncorrected:
correction_type_results_Scaden$uncorrected = Estimate.Proportions('./Data/bulk/CRC/NICs/train_test_sets/train.txt',
                                                                  './Data/CRC_reference/UB_Scaden.h5ad', 'Deconv_cellTypes', 'Scaden',
                                                                  dir.results=tempdir(), datasets='', min.expression=5,
                                                                  batch.size=as.integer(50), learning.rate=0.001, n.steps=as.integer(100))
rownames(correction_type_results_Scaden$uncorrected) = names(cell_types_scaden)
invisible(gc())
# - corrected Before:
correction_type_results_Scaden$corrected.Before = Estimate.Proportions('./Data/bulk/CRC/NICs/train_test_sets/train.txt',
                                                                       './Data/CRC_reference/B_Scaden.h5ad', 'Deconv_cellTypes',
                                                                       'Scaden', dir.results=tempdir(), datasets='',
                                                                       min.expression=0.1, batch.size=as.integer(128), learning.rate=0.001,
                                                                       n.steps=as.integer(500))
rownames(correction_type_results_Scaden$corrected.Before) = names(cell_types_scaden)
invisible(gc())
# - corrected After:
res_temp = Estimate.Proportions('./Data/bulk/CRC/NICs/train_test_sets/train.txt', './Data/CRC_reference/UB_Scaden.h5ad', 'Deconv_cellTypes',
                                'Scaden', dir.results=tempdir(), datasets='', min.expression=0.1, batch.size=as.integer(50),
                                learning.rate=0.001, n.steps=as.integer(1000))
res_temp = correct_fractions_mRNABias(res_temp, rna_correction_vec_scaden)
rownames(res_temp) = names(cell_types_scaden)
correction_type_results_Scaden$corrected.After = res_temp
invisible(gc())
# 3.7.2. Save results:
saveRDS(correction_type_results_Scaden,
        './Results/1_Parameter_Optimization/train/correction_type_analysis_res/Scaden.Rdata')
# 3.7.3. Plot RMSE vs correlation by each cell-types, coloured by correction type:
scatterPlot_rmsecorr_correction_celltype(correction_type_results_Scaden, ground_truth_proportions, train_samples,
                                         ylim=c(-.4,1), xlim=c(0,.45))


# 3.8. SCDC
correction_type_results_SCDC = list()
# Best combo: LAD.uncorrected_A1B1
which.min(train_results$SCDC['MEAN',])
# Best combo for correction before estimation: LAD.corrected.Before_A1B1 
which.min(train_results$SCDC['MEAN',grep('corrected.Before', colnames(train_results$SCDC))])
# Best combo for correction after estimation: LAD.corrected.After_A1B1
which.min(train_results$SCDC['MEAN',grep('corrected.After', colnames(train_results$SCDC))])
# 3.8.1. Run SCDC for LAD.uncorrected_A1B1, LAD.corrected.Before_A1B1 and LAD.corrected.After_A1B1
# - uncorrected:
res_temp = Estimate.Proportions(train, references$uncorrected, 'Deconv_cellTypes', 'SCDC', scdc.method='multiple', n.iterations=100,
                                m.use.gridSearch=FALSE, m.search.length=0.01, m.dataset.var='Dataset', m.weights.method='LAD')[['LAD']]
correction_type_results_SCDC$uncorrected = res_temp
invisible(gc())
# - corrected Before:
correction_type_results_SCDC$corrected.Before = Estimate.Proportions(train, references$corrected, 'Deconv_cellTypes', 'SCDC',
                                                                     scdc.method='multiple', n.iterations=100,
                                                                     m.use.gridSearch=FALSE, m.search.length=0.01,
                                                                     m.dataset.var='Dataset',
                                                                     m.weights.method='LAD')[['LAD']]
invisible(gc())
# - corrected After:
correction_type_results_SCDC$corrected.After = correct_fractions_mRNABias(res_temp, rna_correction_vec)
invisible(gc())
# 3.8.2. Save results:
saveRDS(correction_type_results_SCDC,
        './Results/1_Parameter_Optimization/train/correction_type_analysis_res/SCDC.Rdata')
# 3.8.3. Plot RMSE vs correlation by each cell-types, coloured by correction type:
scatterPlot_rmsecorr_correction_celltype(correction_type_results_SCDC, ground_truth_proportions, train_samples, ylim=c(-.1,1), xlim=c(0,.6))


