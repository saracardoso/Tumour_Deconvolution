
# ---
# - Organize the different combinations of parameters possible
# ---


# Initialize the parameter combinations

parameters_to_optimize=list()

# - CIBERSORTx
parameters_to_optimize$CIBERSORTx = expand.grid(c('YES', 'NO'), c(0, 50, 100, 500, 1000), c(0, 0.5), stringsAsFactors=F)
rownames(parameters_to_optimize$CIBERSORTx) = apply(expand.grid(c('A1', 'A2'), c('B1', 'B2', 'B3', 'B4', 'B5'), c('C1', 'C2'),
                                                                stringsAsFactors=F), 1, paste, collapse='')
colnames(parameters_to_optimize$CIBERSORTx) = c('batch_correction', 'permutations', 'min_expression')

# - AutoGeneS
parameters_to_optimize$AutoGeneS = expand.grid(c('nusvr', 'nnls', 'linear'), c('YES', 'NO'), c('fixed', 'standard'),
                                               c(100, 500, 1000), c(200, 400, 'all'), c('c', 'cd', 'd'),
                                               stringsAsFactors=F)
rownames(parameters_to_optimize$AutoGeneS) = apply(expand.grid(c('A1', 'A2', 'A3'), c('B1', 'B2'), c('C1', 'C2'),
                                                               c('D1', 'D2', 'D3'), c('E1', 'E2', 'E3'), c('F1', 'F2', 'F3'),
                                                               stringsAsFactors=F), 1, paste, collapse='')
colnames(parameters_to_optimize$AutoGeneS) = c('model', 'gene.markers', 'mode', 'n.iterations', 'n.genes', 'weights')
# n.genes is only important if mode is 'fixed'
standard_combos = !duplicated(parameters_to_optimize$AutoGeneS[,-5]) & parameters_to_optimize$AutoGeneS$mode=='standard'
parameters_to_optimize$AutoGeneS = rbind(parameters_to_optimize$AutoGeneS[parameters_to_optimize$AutoGeneS$mode=='fixed',],
                                         parameters_to_optimize$AutoGeneS[standard_combos,])
# When a set of marker genes are given (parameter B1), the number of genes in fixed mode can only go up to ... (nº of markers in our
# reference). Therefore, when B1 and C1, n.genes will be 200, 400 or all genes in markers. If markers are not given 'all' will be 800.

# - MuSiC_woGrouping
parameters_to_optimize$MuSiC_woGrouping = expand.grid(c(100, 500, 1000), c(TRUE, FALSE), c(TRUE, FALSE),
                                           c('weighted', 'all genes'), c('YES', 'NO'),
                                           stringsAsFactors=F)
rownames(parameters_to_optimize$MuSiC_woGrouping) = apply(expand.grid(c('A1', 'A2', 'A3'), c('B1', 'B2'), c('C1', 'C2'),
                                                           c('D1', 'D2'), c('E1', 'E2'),
                                                           stringsAsFactors=F), 1, paste, collapse='')
colnames(parameters_to_optimize$MuSiC_woGrouping) = c('n.iterations', 'center', 'normalise', 'wo.method', 'wo.gene.markers')

# - MuSiC_wGrouping
parameters_to_optimize$MuSiC_wGrouping = expand.grid(c(100, 500, 1000), c(TRUE, FALSE), c(TRUE, FALSE),
                                                     stringsAsFactors=F)
rownames(parameters_to_optimize$MuSiC_wGrouping) = apply(expand.grid(c('A1', 'A2', 'A3'), c('B1', 'B2'), c('C1', 'C2'),
                                                                     stringsAsFactors=F), 1, paste, collapse='')
colnames(parameters_to_optimize$MuSiC_wGrouping) = c('n.iterations', 'center', 'normalise')

# - BSeqSC
parameters_to_optimize$BSeqSC = expand.grid(c('YES', 'NO'), c(0, 50, 100, 500, 1000), stringsAsFactors=F)
rownames(parameters_to_optimize$BSeqSC) = apply(expand.grid(c('A1', 'A2'), c('B1', 'B2', 'B3', 'B4', 'B5'),
                                                            stringsAsFactors=F), 1, paste, collapse='')
colnames(parameters_to_optimize$BSeqSC) = c('ct_scale', 'permutations')

# - SCDC
parameters_to_optimize$SCDC = expand.grid(c(0.01, 0.05, 0.1), c(100, 500, 1000),
                                          stringsAsFactors=F)
rownames(parameters_to_optimize$SCDC) = apply(expand.grid(c('A1', 'A2', 'A3'), c('B1', 'B2', 'B3'),
                                                          stringsAsFactors=F), 1, paste, collapse='')
colnames(parameters_to_optimize$SCDC) = c('m.search.length', 'm.n.iterations')

# - DWLS
parameters_to_optimize$DWLS = expand.grid(c(0.5, 0.8, 1), c(0.01, 0.05),
                                          stringsAsFactors=F)
rownames(parameters_to_optimize$DWLS) = apply(expand.grid(c('A1', 'A2', 'A3'), c('B1', 'B2'),
                                                          stringsAsFactors=F), 1, paste, collapse='')
colnames(parameters_to_optimize$DWLS) = c('diff.cutoff', 'pval.cutoff')

# - BisqueRNA
parameters_to_optimize$BisqueRNA = expand.grid(c('YES', 'NO'), c(TRUE, FALSE),
                                               stringsAsFactors=F)
rownames(parameters_to_optimize$BisqueRNA) = apply(expand.grid(c('A1', 'A2'), c('B1', 'B2'),
                                                               stringsAsFactors=F), 1, paste, collapse='')
colnames(parameters_to_optimize$BisqueRNA) = c('gene.markers', 'old.cpm')

# - MOMF
parameters_to_optimize$MOMF = expand.grid(c('KL', 'IS'), c(1, 2, 5), c(100, 500, 1000),
                                          stringsAsFactors=F)
rownames(parameters_to_optimize$MOMF) = apply(expand.grid(c('A1', 'A2'), c('B1', 'B2', 'B3'), c('C1', 'C2', 'C3'),
                                                          stringsAsFactors=F), 1, paste, collapse='')
colnames(parameters_to_optimize$MOMF) = c('method', 'rho', 'num.iter')

# - Scaden
parameters_to_optimize$Scaden = expand.grid(c(0.1, 1, 5), c(50, 128, 200), c(0.0001, 0.001, 0.01),
                                            c(100, 500, 1000),
                                            stringsAsFactors=F)
rownames(parameters_to_optimize$Scaden) = apply(expand.grid(c('A1', 'A2', 'A3'), c('B1', 'B2', 'B3'),
                                                            c('C1', 'C2', 'C3'), c('D1', 'D2', 'D3'),
                                                            stringsAsFactors=F), 1, paste, collapse='')
colnames(parameters_to_optimize$Scaden) = c('min.expression', 'batch.size', 'learning.rate', 'n.steps')

# - DigitalDLSorter
parameters_to_optimize$DigitalDLSorter = expand.grid(c(0, 1), c('no', 'limit_1000'), c(110, 1100), c('YES', 'NO'), c(100, 500),
                                                     c(10000, 20000), c('YES', 'NO'), c('both', 'single-cell', 'bulk'), c('YES', 'NO'),
                                                     stringsAsFactors=F)
rownames(parameters_to_optimize$DigitalDLSorter) = apply(expand.grid(c('A1', 'A2'), c('B1', 'B2'), c('C1', 'C2'), c('D1', 'D2'),
                                                                     c('E1', 'E2'), c('F1', 'F2'), c('G1', 'G2'),
                                                                     c('H1', 'H2', 'H3'), c('I1', 'I2'),
                                                                     stringsAsFactors=F), 1, paste, collapse='')
colnames(parameters_to_optimize$DigitalDLSorter) = c('load_min.counts', 'simulation', 'simul_subset.cells',
                                                     'simul_proportional', 'bulk_n.cells', 'bulk_num.bulk.samples',
                                                     'bulk_balanced.type.cells', 'train_combine', 'predict_normalize')

# Save combinations in csv files
for(method in names(parameters_to_optimize)){
  write.csv(parameters_to_optimize[[method]],
            paste('./Results/1_Parameter_Optimization/parameters/', method, '.csv', sep=''))
}

