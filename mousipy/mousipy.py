import pandas as pd
import numpy as np
import scanpy as sc
from tqdm import tqdm

h2m_tab = pd.read_csv('mousipy/biomart/human_to_mouse_biomart_export.csv').set_index('Gene name')
m2h_tab = pd.read_csv('mousipy/biomart/mouse_to_human_biomart_export.csv').set_index('Gene name')

def check_orthologs(var_names, tab=m2h_tab):
    direct = {}
    multiple = {}
    no_hit = []
    no_index = []
    for gene in tqdm(var_names):
        if gene in tab.index:
            # gene in index
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
    # direct hits can be used as is
    # for those with no entry in the database we assume that uppercase-ing works out as a good guess
    # (after excluding some gene symbols which usually do not map)
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

def translate_multiple(adata, original_data, multiple):
    # multiple hits are tricky:
    # we will add the counts of a mouse gene to ALL the homologs it maps to.
    # Unfortunately this takes pretty long (~10mins for ~1000 genes)
    from scipy.sparse import csr_matrix, hstack
    X = adata.X
    var = adata.var
    for mgene, hgenes in tqdm(multiple.items()):
        for hgene in hgenes:
            if hgene not in list(var.index):
                # Add counts to new gene
                X = csr_matrix(hstack((X, original_data[:, mgene].X)))
                var.loc[hgene] = None
                var.loc[hgene, 'original_gene_symbol'] = 'multiple'
            else:
                # Add counts to existing gene
                idx = np.where(np.array(list(var.index))==hgene)[0][0]
                X[:, idx] += original_data[:, mgene].X
    return sc.AnnData(X, adata.obs, var, adata.uns, adata.obsm)

def translate(adata, tab=m2h_tab):
    direct, multiple, no_hit, no_index = check_orthologs(adata.var_names, tab=m2h_tab)

    # for those with an entry but no ortholog we assume that there really is no known ortholog
    # which means we ignore genes in m2h_no_hit
    bdata = translate_direct(adata, direct, no_index)
    bdata = translate_multiple(bdata, adata, multiple)
    return bdata

def test():
    import scvelo
    adata = scv.datasets.pancreas()
    humanized_pancreas = translate(adata)
