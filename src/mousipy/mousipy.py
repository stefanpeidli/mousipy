import os
from re import search

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix, issparse
from tqdm import tqdm

# Biomart tables
path = os.path.abspath(os.path.dirname(__file__))
h2m_tab_biomart = pd.read_csv(os.path.join(path, './biomart/human_to_mouse.csv')).set_index('Gene name')
m2h_tab_biomart = pd.read_csv(os.path.join(path, './biomart/mouse_to_human.csv')).set_index('Gene name')

# HCOP tables
h2m_tab_hcop = pd.read_csv(os.path.join(path, './hcop/human_to_mouse.csv')).set_index('Gene name')
m2h_tab_hcop = pd.read_csv(os.path.join(path, './hcop/mouse_to_human.csv')).set_index('Gene name')


def make_dense(X):
    # robustly make an array dense
    if issparse(X):
        return X.A
    else:
        return X


def identify_format_and_organism(txt):
    '''Identify the format and organism of a gene symbol.
    Currently only supports mouse and human genes.
    Parameters
    ----------
    txt : str
        A gene string.

    Returns
    -------
    format: str
        Either 'Ensembl' or 'Symbol'.
    organism: str
        Either 'Mouse' or 'Human'.
    '''
    if search('ENSG[0-9]{11}', txt):
        format = 'Ensembl'
        organism = 'Human'
    elif search('ENSMUSG[0-9]{11}', txt):
        format = 'Ensembl'
        organism = 'Mouse'
    else:
        format = 'Symbol'
        if txt.isupper():  # probably not very robust
            organism = 'Human'
        else:
            organism = 'Mouse'
    return format, organism


def check_orthologs(var_names, target=None, tab=None, verbose=False, source='biomart'):
    """Check for orthologs from a list of gene symbols in a biomart table.
    In essence, this function goes through the list of input genes and checks
    if they have a single ortholog in the table. If they have multiple orthologs,
    no ortholog, or are not found in the table, they are added to the respective
    return dictionary or list.
    
    Parameters
    ----------
    var_names : list or list-like
        A list of (mouse) gene symbols.
    target : str (either 'Ensembl' or 'Symbol') or None, optional
        The target gene names that will be translated to. 
        If None, translates Ensembl to Ensembl and Symbol to Symbol.
    tab :
        A biomart table. If None, the default mouse-to-human table is used.
    verbose : bool, optional
        Whether to print the number progress of the translation. The default is False.

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
    if target not in ['Ensembl', 'Symbol', None]:
        raise ValueError('target must be either "Ensembl", "Symbol", or None')

    source_lower = source.lower()

    if source_lower == 'hcop':
        h2m_tab = h2m_tab_hcop
        m2h_tab = m2h_tab_hcop
    elif source_lower == 'biomart':
        h2m_tab = h2m_tab_biomart
        m2h_tab = m2h_tab_biomart
    else:
        raise ValueError('source must be either "HCOP" or "biomart"')

    input_format, input_organism = identify_format_and_organism(var_names[0])
    output_organism = 'Human' if input_organism == 'Mouse' else 'Mouse'
    if not isinstance(tab, pd.DataFrame):
        tab = m2h_tab if input_organism == 'Mouse' else h2m_tab
    if input_format == 'Ensembl':
        tab = tab.reset_index().set_index('Gene stable ID')
        target = 'Ensembl' if target is None else target
    else:
        # tab is already set to gene symbols as index
        target = 'Symbol' if target is None else target
    target_key = f'{output_organism} gene stable ID' if target == 'Ensembl' else f'{output_organism} gene name'

    # very slow, might want to speed up. Maybe with a lookup table?
    direct = {}  # genes with a single ortholog
    multiple = {}  # genes with multiple orthologs
    no_hit = []  # genes with no ortholog
    no_index = []  # genes not in the index of the table
    fct = tqdm if verbose else lambda x: x
    for gene in fct(var_names):
        if gene in tab.index:
            # gene found in index of the table
            x = tab[target_key].loc[gene]
            if isinstance(x, pd.Series):
                # multiple hits
                vals = pd.unique(x.values)
                vals = vals[~pd.isna(vals)]
                if len(vals) > 1:
                    # multiple actual hits
                    multiple[gene] = vals
                elif len(vals) == 1:
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


def translate_multiple(adata, original_data, multiple, stay_sparse=False, verbose=False):
    """
    Adds the counts of multiple-hit genes to ALL their orthologs.
    """
    var = adata.var.copy()
    ortholog_indices = {gene: i for i, gene in enumerate(var.index)}

    if stay_sparse:
        # Sparse implementation remains unchanged
        X = adata.X.copy()
        for mgene, hgenes in (tqdm(multiple.items()) if verbose else multiple.items()):
            mgene_data = make_dense(original_data[:, mgene].X)

            for hgene in hgenes:
                if hgene not in ortholog_indices:
                    # Create a new DataFrame row for the new gene
                    new_row = pd.DataFrame({col: pd.NA for col in var.columns}, index=[hgene])
                    new_row['original_gene_symbol'] = 'multiple'
                    var = pd.concat([var, new_row])

                    X = csr_matrix(np.hstack((X.toarray(), mgene_data.reshape(-1, 1))))
                    ortholog_indices[hgene] = X.shape[1] - 1
                else:
                    idx = ortholog_indices[hgene]
                    X[:, idx] += csr_matrix(mgene_data).reshape(-1, 1)
    else:
        # Dense implementation
        num_new_genes = sum(1 for hgenes in multiple.values() for hgene in hgenes if hgene not in ortholog_indices)
        X = make_dense(adata.X)
        new_data = np.zeros((X.shape[0], X.shape[1] + num_new_genes))

        new_data[:, :X.shape[1]] = X
        next_new_gene_idx = X.shape[1]

        for mgene, hgenes in (tqdm(multiple.items()) if verbose else multiple.items()):
            mgene_data = make_dense(original_data[:, mgene].X).reshape(-1, 1)

            for hgene in hgenes:
                if hgene not in ortholog_indices:
                    # Create a new DataFrame row for the new gene
                    new_row = pd.DataFrame({col: pd.NA for col in var.columns}, index=[hgene])
                    new_row['original_gene_symbol'] = 'multiple'
                    var = pd.concat([var, new_row])

                    new_data[:, next_new_gene_idx] = mgene_data.ravel()
                    ortholog_indices[hgene] = next_new_gene_idx
                    next_new_gene_idx += 1
                else:
                    idx = ortholog_indices[hgene]
                    new_data[:, idx] += mgene_data.ravel()

        X = new_data

    # Check the dimensions of X and var
    if X.shape[1] != var.shape[0]:
        # If they do not match, modify var to match the dimensions
        missing_rows = X.shape[1] - var.shape[0]
        additional_rows = pd.DataFrame(index=range(var.shape[0], X.shape[1]))
        var = pd.concat([var, additional_rows])

    return AnnData(X, adata.obs, var, adata.uns, adata.obsm)


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


def translate(adata, target=None, stay_sparse=False, verbose=True, source='biomart'):
    """Translates adata.var between mouse and human using orthologs from biomart.
    Accepts either Ensembl or Symbol gene names as adata.var_names and can translate
    to either Ensembl or Symbol gene names. If the target is None, translates
    to the same gene name type as the input.
    
    Parameters
    ----------
    adata : AnnData object
        Single cell object to translate.
    target : str (either 'Ensembl' or 'Symbol') or None, optional
        The target gene names that will be translated to. 
        If None, translates Ensembl to Ensembl and Symbol to Symbol.
    stay_sparse :
        If True, forbids this function to locally densify the count matrix.

    Returns
    -------
    AnnData object
        The adata with translated and unique features.
    """
    source_lower = source.lower()

    if source_lower not in ['hcop', 'biomart']:
        raise ValueError('source must be either "HCOP" or "biomart"')

    direct, multiple, no_hit, no_index = check_orthologs(adata.var_names, tab=None, target=target, verbose=verbose, source=source_lower)
    if verbose:
        print('Found direct orthologs for {} genes.'.format(len(direct)))
        print('Found multiple orthologs for {} genes.'.format(len(multiple)))
        print('Found no orthologs for {} genes.'.format(len(no_hit)))
        print('Found no index in biomart for {} genes.'.format(len(no_index)))

    # for those with an entry but no ortholog we assume that there really is no known ortholog
    # which means we ignore genes in m2h_no_hit
    bdata = translate_direct(adata, direct, no_index)
    bdata = translate_multiple(bdata, adata, multiple, stay_sparse=stay_sparse, verbose=verbose)
    bdata = collapse_duplicate_genes(bdata, stay_sparse=stay_sparse)
    return bdata
