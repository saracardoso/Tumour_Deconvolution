#########################################################################
## Simulation of artificial bulk RNA-seq datasets from scRNA-seq data   #
########################################################################

import pandas as pd
import numpy as np
import glob
import os
import argparse
import gc

def create_fractions(no_celltypes):
    """
    Create random fractions
    :param no_celltypes: number of fractions to create
    :return: list of random fracs of length no_cellttypes
    """
    fracs = np.random.rand(no_celltypes)
    fracs_sum = np.sum(fracs)
    fracs = np.divide(fracs, fracs_sum)
    return fracs


def create_subsample(x, y, sample_size, celltypes, available_celltypes, sparse=False):
    """
    Generate artifical bulk subsample with random fractions of celltypes
    If sparse is set to true, add random celltypes to the missing celltypes
    :param cells: scRNA-seq expression matrix
    :param labels: celltype labels of single cell data (same order)
    :param sample_size: number of cells to use
    :param celltypes: all celltypes available in all datasets
    :param available_celltypes: celltypes available in currently processed dataset
    :param sparse: whether to create a sparse datasets (missing celltypes)
    :return: expression of artificial subsample along with celltype fractions
    """

    if sparse:
        no_keep = np.random.randint(1, len(available_celltypes))
        keep = np.random.choice(list(range(len(available_celltypes))), size=no_keep, replace=False)
        available_celltypes = [available_celltypes[i] for i in keep]

    no_avail_cts = len(available_celltypes)

    # Create fractions for available celltypes
    fracs = create_fractions(no_celltypes=no_avail_cts)
    samp_fracs = np.multiply(fracs, sample_size)
    samp_fracs = list(map(int, samp_fracs))

    # Make complete fracions
    fracs_complete = [0] * len(celltypes)
    for i,act in enumerate(available_celltypes):
        idx = celltypes.index(act)
        fracs_complete[idx] = fracs[i]

    artificial_samples = []
    for i in range(no_avail_cts):
        ct = available_celltypes[i]
        cells_sub = x.loc[np.array(y['Celltype'] == ct),:]
        cells_fraction = np.random.randint(0, cells_sub.shape[0], samp_fracs[i])
        cells_sub = cells_sub.iloc[cells_fraction, :]
        artificial_samples.append(cells_sub)

    df_samp = pd.concat(artificial_samples, axis=0)
    df_samp = df_samp.sum(axis=0)

    return (df_samp, fracs_complete)


def create_subsample_dataset(x, y, sample_size, celltypes, no_samples):
    """
    Generate many artifial bulk samples with known fractions
    This function will create normal and sparse samples (no_samples)
    :param cells:
    :param labels:
    :param sample_size:
    :param no_celltypes:
    :param no_samples:
    :return: dataset and corresponding labels
    """
    X = []
    Y = []
    
    available_celltypes = list(set(y['Celltype'].tolist()))

    # Create normal samples
    #print('Creating normal samples...')
    for i in range(no_samples):
        sample, label = create_subsample(x, y, sample_size, celltypes, available_celltypes)
        X.append(sample)
        Y.append(label)

    # Create sparse samples
    n_sparse = int(no_samples)
    #n_sparse = 0
    #print('Creating sparse samples...')
    for i in range(n_sparse):
        sample, label = create_subsample(x, y, sample_size, celltypes, available_celltypes, sparse=True)
        X.append(sample)
        Y.append(label)
    X = pd.concat(X, axis=1).T
    Y = pd.DataFrame(Y, columns=celltypes)

    # Shuffle
    #X, Y = shuffle_dataset(X, Y)

    return (X, Y)

def filter_for_celltypes(x, y, celltypes):
    """
    Filter data for cells belonging to specified celltypes
    :param x:
    :param y:
    :param celltypes:
    :return:
    """
    cts = list(y['Celltype'])
    keep = [elem in celltypes for elem in cts]
    x = x.loc[keep, :]
    y = y.loc[keep, :]
    return (x, y)



def shuffle_dataset(x, y):
    """
    Shuffle dataset while keeping x and y in synch
    :param x:
    :param y:
    :return:
    """
    idx = np.random.permutation(x.index)
    x_shuff = x.reindex(idx)
    y_shuff = y.reindex(idx)
    return (x_shuff, y_shuff)

def filter_matrix_signature(mat, genes):
    """
    Filter expression matrix using given genes
    Accounts for the possibility that some genes might not be in the matrix, for these genes
    a column with zeros is added to the matrix
    :param mat:
    :param genes:
    :return: filtered matrix
    """
    n_cells = mat.shape[0]
    avail_genes = mat.columns
    filt_genes = [g for g in genes if g in avail_genes]
    missing_genes = [g for g in genes if g not in avail_genes]
    mat = mat[filt_genes]
    for mg in missing_genes:
        mat[mg] = np.zeros(n_cells)
    mat = mat[genes]
    return mat

def load_dataset(name, dir, pattern):
    """
    Load a dataset given its name and the directory
    :param name: name of the dataset
    :param dir: directory containing the data
    :param sig_genes: the signature genes for filtering
    :return: X, Y
    """
    pattern = pattern.replace("*", "")
    print("Loading " + name + " dataset ...")
    x = pd.read_table(dir + name + pattern, index_col=0)
    y = pd.read_table(dir + name + "_celltypes.txt")
    return (x, y)

def merge_unkown_celltypes(y, unknown_celltypes):
    """
    Merge all unknown celltypes together
    :param x:
    :param y:
    :param unknown_celltypes:
    :return:
    """
    celltypes = list(y['Celltype'])
    new_celltypes = ["Unknown" if x in unknown_celltypes else x for x in celltypes]
    y['Celltype'] = new_celltypes
    return y

def collect_celltypes(ys):
    """
    Collect all available celltypes given all dataset labels
    :param ys: list of dataset labels
    :return: list of available celltypes
    """
    ct_list = [list(set(y['Celltype'].tolist())) for y in ys]
    celltypes = set()
    for ct in ct_list:
        celltypes = celltypes.union(set(ct))
    celltypes = list(celltypes)
    return celltypes

def get_common_genes(xs, type='intersection'):
    """
    Get common genes for all matrices xs
    Can either be the union or the intersection (default) of all genes
    :param xs: cell x gene matrices
    :return: list of common genes
    """
    genes = []
    for x in xs:
        genes.append(list(x.columns))

    genes = [set(g) for g in genes]
    com_genes = genes[0]
    if type=='union':
        for gi in range(1, len(genes)):
            com_genes = com_genes.union(genes[gi])
    elif type=='intersection':
        for gi in range(1, len(genes)):
            com_genes = com_genes.intersection(genes[gi])

    else:
        exit('Wrong type selected to get common genes. Exiting.')

    if len(com_genes) == 0:
        exit("No common genes found. Exiting.")

    return list(com_genes)

def generate_signature(x, y):
    """
    Generate signature of matrix using celltypes y
    :param x: expression matrix
    :param y: celltypes
    :return: mean-expression per celltype
    """

    signature_matrix = []
    celltypes = list(set(y['Celltype']))
    for ct in celltypes:
        ct_exp = x.loc[np.array(y['Celltype'] == ct), :]
        ct_exp = ct_exp.mean(axis=0)
        signature_matrix.append(ct_exp)

    signature_matrix = pd.concat(signature_matrix, axis=1)
    signature_matrix.columns = celltypes
    return signature_matrix
  
