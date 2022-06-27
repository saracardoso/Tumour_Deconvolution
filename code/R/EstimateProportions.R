algorithms_file = './code/R/EstimateProportions_Algorithms.R'
source(algorithms_file)

#---
#- MAIN DECONVOLUTION FUNCTION FOR ESTIMATING PROPORTIONS
#---

# Methods available: MuSiC_woGrouping, MuSiC_wGrouping, SCDC_single, SCDC_multiple, BisqueRNA, DWLS, AutoGeneS, MOMF, CPM, Scaden
# Note - CIBERSORTx: web tool, only plot functions can be used to visualize results. If one wants to use more than one single-cell reference,
#                    use scRefs_batch_correction to generate the 'raw' matrix to submit to CIBERSORTx.
# Note - Bseq-SC: uses CIBERSORT to calculate proportions. Thus, reference matrix creation is available through create_reference function.
#                 To calculate proportions, one should submit the reference generated and the bulk data to CIBERSORTx. Then, plot functions
#                 can be used to visualize the results.
Estimate.Proportions = function(bulk.data, sc.reference, cellType.var, method, ...){
  methods = c('DWLS', 'MOMF', 'AutoGeneS', 'MuSiC', 'SCDC', 'BisqueRNA', 'Scaden', 'DigitalDLSorter')
  if(!method%in%methods) stop(paste('Invalid method! Available methods:', methods, collapse='; '))
  class(bulk.data) = method
  proportions = estimate_proportions(bulk.data, sc.reference, cellType.var, ...)
  return(proportions)
}



estimate_proportions = function(bulk.data, sc.reference, cellType.var, ...)
  UseMethod('estimate_proportions')

estimate_proportions.DWLS = function(bulk.data, sc.reference, cellType.var, diff.cutoff=0.5, pval.cutoff=0.01, rna.bias=FALSE, bias.vec=NULL){
  proportions = deconv_dwls(bulk.data$data, sc.reference, cellType.var, diff.cutoff, pval.cutoff, rna.bias, bias.vec)
  return(proportions)
}

estimate_proportions.MOMF = function(bulk.data, sc.reference, cellType.var, momf.method='KL', rho=2, num.iter=5000, rna.bias=FALSE, bias_vec=NULL){
  proportions = deconv_momf(bulk.data$data, sc.reference, cellType.var, momf.method, rho, num.iter, rna.bias, bias.vec)
  return(proportions)
}

estimate_proportions.AutoGeneS = function(bulk.data, sc.reference, cellType.var, model='nnls', gene.markers=NULL, mode='fixed',
                                          n.iterations=5000L, n.genes=400L, weights=list(-1,1)){
  proportions = deconv_autogenes(bulk.data$data, sc.reference, cellType.var, model, gene.markers, mode, n.iterations, n.genes, weights)
  return(proportions)
}

estimate_proportions.MuSiC = function(bulk.data, sc.reference, cellType.var, music.method='without.Grouping', n.iterations=1000, center=FALSE,
                                      normalize=FALSE,
                                      wo.method='Est.prop.weighted', wo.gene.markers=NULL, wo.bias.vec=NULL,
                                      w.celltypes.Cluster.List=NULL, w.markers.log2FC=1,
                                      w.markers.pval=0.01, w.markers.only.Pos=FALSE
                                      ){
  scReference = prepare_expressionSet_scReferences(sc.reference, cellType.var)
  bulk_data_expressionSet = Biobase::ExpressionSet(data.matrix(bulk.data$data))
  if(music.method=='without.Grouping')
    proportions = deconv_music_woGrouping(bulk_data_expressionSet, scReference, wo.method, wo.gene.markers,
                                          n.iterations, center, normalize, wo.bias.vec)
  else if(music.method=='with.Grouping')
    proportions = deconv_music_wGrouping(bulk_data_expressionSet, scReference, w.celltypes.Cluster.List, w.markers.log2FC,
                                         w.markers.pval, w.markers.only.Pos, n.iterations, center, normalize)
  
  return(proportions)
}

estimate_proportions.SCDC = function(bulk.data, sc.reference, cellType.var, scdc.method='multiple', n.iterations=1000,
                                      m.use.gridSearch=FALSE, m.search.length=0.05, m.dataset.var='', m.weights.method=NULL){
  bulk_data_expressionSet = Biobase::ExpressionSet(data.matrix(bulk.data$data))
  if(scdc.method=='single'){
    scReference = prepare_expressionSet_scReferences(sc.reference, cellType.var)
    proportions = deconv_scdc_single_ref(bulk_data_expressionSet, scReference, n.iterations)
  }
  else if(scdc.method=='multiple'){
    if(!m.dataset.var%in%colnames(sc.reference$metadata)) stop('Invalid m.dataset.var.')
    cat('Creating multiple scReference datasets acording to variable', m.dataset.var, '\n')
    scReference = list()
    for(var_dat in unique(sc.reference$metadata[,m.dataset.var])){
      var_dat_samples = rownames(sc.reference$metadata)[sc.reference$metadata[,m.dataset.var]==var_dat]
      temp = sc.reference
      temp$data = temp$data[, var_dat_samples]
      temp$metadata = temp$metadata[var_dat_samples, ]
      scReference[[as.character(var_dat)]] = prepare_expressionSet_scReferences(temp, cellType.var)
    }
    proportions = deconv_scdc_multiple_refs(bulk_data_expressionSet, scReference, n.iterations,
                                            m.use.gridSearch, m.search.length, m.weights.method)
  }
  return(proportions)
}

estimate_proportions.BisqueRNA = function(bulk.data, sc.reference, cellType.var, gene.markers=NULL, use.overlap=FALSE, old.cpm=TRUE){
  scReference = prepare_expressionSet_scReferences(sc.reference, cellType.var)
  bulk_data_expressionSet = Biobase::ExpressionSet(data.matrix(bulk.data$data))
  proportions = deconv_bisquerna(bulk_data_expressionSet, scReference, gene.markers, use.overlap, old.cpm)
  return(proportions)
}

estimate_proportions.Scaden = function(bulk.data, sc.reference, cellType.var, dir.results=tempdir(), datasets='', min.expression=0.1,
                                       batch.size=128L, learning.rate=0.0001, n.steps=1000L){
  proportions = deconv_scaden(bulk.data, sc.reference, cellType.var, dir.results, datasets, min.expression, batch.size, learning.rate, n.steps)
  return(proportions)
}

estimate_proportions.DigitalDLSorter = function(bulk.data, sc.reference, cellType.var, sample.var = 'Sample',
                                                pipeline='train_predict', simulate.cells=TRUE, return.model=FALSE,
                                                load_min.counts=0, load_min.cells=0,
                                                simul_subset.cells=NULL, simul_proportional=TRUE, simul_n.cells=100, simul_cell.types=NULL,
                                                simul_add.method='fixed',
                                                bulk_n.cells=100, bulk_balanced.type.cells=FALSE, bulk_num.bulk.samples=250,
                                                train_combine='both', train_num.epocs=10, train_view.metrics.plot=FALSE,
                                                predict_normalize=TRUE, batch.size=64, threads=1){
  if(!pipeline%in%c('train', 'predict', 'train_predict')) stop('Invalid pipeline! Available pipelines: train, predict, train_predict')
  if(pipeline%in%c('train', 'train_predict')){
    message('Trainning...')
    DDLSorter_obj = prepare_singlecellexperiment_scReferences(sc.reference, cellType.var, load_min.cells, load_min.counts, 'DDLS')
    # In this case, sc.reference is the original reference data!
    if(simulate.cells){
      DDLSorter_obj = simulate_cells(DDLSorter_obj, simul_n.cells, simul_cell.types, simul_subset.cells, simul_proportional, threads,
                                     simul_add.method)
      invisible(gc())
    }
    # Calculate minimum and maximum proportions of cell-types by sample in reference to create the minmaxprobs argument
    cts = unique(sc.reference$metadata[,cellType.var])
    min_props = rep(1, length(cts))
    max_props = rep(0, length(cts))
    samps = unique(sc.reference$metadata[,sample.var])
    for(samp in samps){
      tab_cts = table(sc.reference$metadata[,cellType.var][sc.reference$metadata[,sample.var]==samp])
      samp_cts = rep(0, length(cts))
      names(samp_cts) = cts
      samp_cts[names(tab_cts)] = tab_cts / sum(tab_cts)
      min_props[samp_cts < min_props] = samp_cts[samp_cts < min_props]
      max_props[samp_cts > max_props] = samp_cts[samp_cts > max_props]
    }
    min_props[min_props==0] = 0.01
    minmaxprobs = data.frame(Cell_Type=cts, from=min_props*100, to=max_props*100)
    # Generate bulk:
    invisible(gc())
    DDLSorter_obj = generate_bulk(DDLSorter_obj, minmaxprobs, bulk_num.bulk.samples, bulk_n.cells, bulk_balanced.type.cells, threads)
    invisible(gc())
    DDLSorter_obj = train_DDLSorter_model(DDLSorter_obj, train_combine, batch.size, train_num.epocs, threads, train_view.metrics.plot)
    invisible(gc())
    if(pipeline == 'train_predict'){
      message('Predicting...')
      proportions = deconv_DigitalDLSorter(DDLSorter_obj, bulk.data, batch.size, predict_normalize, 'deconv')
      invisible(gc())
    }
  }
  else{
    # In this case, sc.reference must be a DigitalDLSorter object!
    DDLSorter_obj = sc.reference
    proportions = deconv_DigitalDLSorter(DDLSorter_obj, bulk.data, batch.size, predict_normalize, 'deconv')
    invisible(gc())
  }
  if(return.model) res_return = list(proportions=proportions, model_DDLS=DDLSorter_obj)
  else res_return = list(proportions=proportions)
  return(res_return)
}
