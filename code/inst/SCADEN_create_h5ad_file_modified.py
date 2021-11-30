"""
Combine artificial bulk datasets and optionally other datasets into h5ad files
usable for scaden training

When using additional datasets, they should be in similar format and best have the same output cell types.
"""

import argparse
import anndata
import glob
import os
import pandas as pd
import numpy as np

"""
Functions
"""

def parse_data(y):
    """
    Parse data and labels and divide them into training and testset
    :param y:
    :return: training and test data and labels
    """
    
    labels = list(y.columns)

    # Transform Y to numpy array and split in train and testset
    yseries = []

    for i in range(y.shape[0]):
        yseries.append(list(y.iloc[i]))
    y = np.array(yseries)

    return y, labels

def sort_celltypes(ratios, labels, ref_labels):
    """
    Bring ratios in correct order of cell types
    :param ratios:
    :param labels:
    :param ref_labels:
    :return:
    """
    idx = [labels.index(x) for x in ref_labels]
    ratios = ratios[:, idx]
    return ratios

def create_AnnData_dataset(xs, ys, celltypes, unknown_cell_types, output_filename):
    
    final_adata = []
    for ds in xs.keys():
      x = xs[ds]
      y = ys[ds]
      
      y, labels = parse_data(y)
      # sort y
      y = sort_celltypes(y, labels, celltypes)
      test = [labels.index(x) for x in celltypes]
      labels = [labels[i] for i in test]
  
      x = x.sort_index(axis=1)
      ratios = pd.DataFrame(y, columns=celltypes)
      ratios['ds'] = pd.Series(np.repeat(ds, y.shape[0]),
                               index=ratios.index)

      ds_anndata = anndata.AnnData(X=x.to_numpy(), obs=ratios, var=pd.DataFrame(columns=[], index=list(x)))
      
      if len(final_adata) == 0:
        final_adata = ds_anndata
      else:
        final_adata.concatenate(ds_anndata)
    
    # add cell types and signature genes
    final_adata.uns['cell_types'] = celltypes
    final_adata.uns['unknown'] = unknown_cell_types
    
    # save data
    final_adata.write(output_filename)
    
    
    
