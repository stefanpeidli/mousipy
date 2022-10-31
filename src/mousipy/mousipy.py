import os
import pandas as pd
import numpy as np
import scanpy as sc

from tqdm import tqdm
from scipy.sparse import csr_matrix, hstack, issparse

# Biomart tables
path = os.path.abspath(os.path.dirname(__file__))
h2m_tab = pd.read_csv(os.path.join(path, './biomart/human_to_mouse_biomart_export.csv')).set_index('Gene name')
m2h_tab = pd.read_csv(os.path.join(path, './biomart/mouse_to_human_biomart_export.csv')).set_index('Gene name')

def make_dense(X):
    # robustly make an array dense
    if issparse(X):
        return X.A
    else:
        return X

def check_orthologs(var_names, tab=None):
    """Check for orthologs from a list of gene symbols in a biomart table.
    Parameters
    ----------
    var_names : list or list-like
        A list of (mouse) gene symbols.
    tab :
        If True, forbids this function to locally densify the count matrix.

    Returns
    -------
    dict
        Dictionary with those input genes mapping to a single ortholog
    dict
        Dictionary of those input genes mapping to multiple orthologs
    list
        List of those input genes found in the table but with no known ortholog
    list
        List of those input genes not found in the table
    """
    tab = tab if isinstance(tab, pd.DataFrame) else m2h_tab

    direct = {}
    multiple = {}
    no_hit = []
    no_index = []
    for gene in tqdm(var_names, leave=False):
        if gene in tab.index:
            # gene found in index of the table
            x = tab['Human gene name'].loc[gene]
            if isinstance(x, pd.Series):
                # multiple hits
                vals = pd.unique(x.values)
                vals = vals[~pd.isna(vals)]
                if len(vals)>1:
                    # multiple actual hits
                    multiple[gene]=vals
                elif len(vals)==1:
                    # one actual hit
                    direct[gene] = vals[0]
                else:
                    # actually no hit (t'was multiple nans)
                    no_hit.append(gene)
            elif pd.isna(x):
                # no hit
                no_hit.append(gene)
            else:
                # one hit
                direct[gene] = x
        else:
            # gene not in index
            no_index.append(gene)
    return direct, multiple, no_hit, no_index

def translate_direct(adata, direct, no_index):
    """Translate all direct hit genes into their orthologs.
    Genes not found in the index of the table will be upper-cased, after
    excluding some gene symbols that usually do not have an ortholog, i.e. genes
    - Starting with 'Gm'
    - Starting with 'RP'
    - Ending with 'Rik'
    - Containing a 'Hist'
    - Containing a 'Olfr'
    - Containing a '.'
    Parameters
    ----------
    adata : AnnData
        AnnData object to translate genes in.
    direct : dict
        Dictionary with those adata genes mapping to a single ortholog
    no_index : list
        List of those adata genes not found in the table.

    Returns
    -------
    AnnData
        Updated original adata.
    """
    guess_genes = [
    x for x in no_index if
    x[:2] != 'Gm' and
    'Rik' not in x and
    x[:2] != 'RP' and
    'Hist' not in x and
    'Olfr' not in x and
    '.' not in x]
    ndata = adata[:, list(direct.keys()) + guess_genes].copy()
    ndata.var['original_gene_symbol'] = list(direct.keys()) + guess_genes
    ndata.var_names = list(direct.values()) + [m.upper() for m in guess_genes]
    return ndata

def translate_multiple(adata, original_data, multiple, stay_sparse=False):
    """Adds the counts of multiple-hit genes to ALL their orthologs.
    Parameters
    ----------
    adata : AnnData
        AnnData object to translate genes in.
    original_data : AnnData
        Original AnnData (before translate_direct).
    multiple : dict
        Dictionary of those adata genes mapping to many orthologs.

    Returns
    -------
    AnnData
        Updated original adata.
    """
    X = adata.X.copy() if stay_sparse else make_dense(adata.X).copy()
    var = adata.var.copy()
    for mgene, hgenes in tqdm(multiple.items(), leave=False):
        for hgene in hgenes:
            if hgene not in list(var.index):
                # Add counts to new gene
                X = np.hstack((X, make_dense(original_data[:, mgene].X)))
                var.loc[hgene] = None
                var.loc[hgene, 'original_gene_symbol'] = 'multiple'
            else:
                # Add counts to existing gene
                idx = np.where(np.array(list(var.index))==hgene)[0][0]
                X[:, [idx]] += make_dense(original_data[:, mgene].X)
    X = X if stay_sparse or not issparse(adata.X) else csr_matrix(X)
    return sc.AnnData(X, adata.obs, var, adata.uns, adata.obsm)

def collapse_duplicate_genes(adata, stay_sparse=False):
    """Collapse duplicate genes by summing up counts to unique entries.
    Adds the counts of duplicate genes to the first duplicate entry.
    Then only keeps the first such entries, making the genes unique without
    changing the number of counts per cell.

    Parameters
    ----------
    adata : AnnData object
        Single cell object to collapse var_names of.
    stay_sparse :
        If True, forbids this function to locally densify the count matrix.

    Returns
    -------
    AnnData object
        The adata with collapsed duplicate genes (has less features now).
    """
    index = adata.var.index
    is_duplicated = index.duplicated()
    duplicated_genes = np.sort(pd.unique(list(index[is_duplicated])))

    if len(duplicated_genes) == 0:
        print('No duplicate genes found. Stopping...')
        return adata

    idxs_to_remove = []
    X = adata.X if stay_sparse else make_dense(adata.X)

    for gene in tqdm(duplicated_genes, leave=False):
        idxs = np.where(index == gene)[0]
        # add later entry counts to first entry counts
        X[:, idxs[0]] += np.sum(X[:, idxs[1:]], axis=1)
        # mark later entries for deletion
        idxs_to_remove += list(idxs[1:])

    adata.X = X if stay_sparse or not issparse(adata.X) else csr_matrix(X)
    return adata[:, np.delete(np.arange(adata.n_vars), idxs_to_remove)].copy()

def translate(adata, stay_sparse=False):
    """Translates adata.var from mouse to human gene symbols using orthologs from biomart.
    Parameters
    ----------
    adata : AnnData object
        Single cell object to translate.
    stay_sparse :
        If True, forbids this function to locally densify the count matrix.

    Returns
    -------
    AnnData object
        The adata with translated and unique features.
    """
    direct, multiple, no_hit, no_index = check_orthologs(adata.var_names, tab=m2h_tab)

    # for those with an entry but no ortholog we assume that there really is no known ortholog
    # which means we ignore genes in m2h_no_hit
    bdata = translate_direct(adata, direct, no_index)
    bdata = translate_multiple(bdata, adata, multiple, stay_sparse=stay_sparse)
    bdata = collapse_duplicate_genes(bdata, stay_sparse=stay_sparse)
    return bdata
