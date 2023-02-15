from mousipy import translate
from scanpy import read

def test_pancreas():
    # download mouse scRNA-seq dataset
    url_datadir = "https://github.com/theislab/scvelo_notebooks/raw/master/"
    url = f"{url_datadir}data/Pancreas/endocrinogenesis_day15.h5ad"
    adata = read("data/Pancreas/endocrinogenesis_day15.h5ad", backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()
    
    humanized_adata = translate(adata)
    assert humanized_adata.n_obs == adata.n_obs, "We lost cells during mapping, which should not happen!"
    assert humanized_adata.n_vars > 10000, "Very few genes (less than 10k) could be mapped! Expecting around 17.5k."

def test_human_PBMC():
    # download human scRNA-seq dataset
    url = 'http://falexwolf.de/data/pbmc3k_raw.h5ad'
    adata = read("data/Pancreas/pbmc3k_raw.h5ad", backup_url=url)
    adata.var_names_make_unique()
    
    mousified_adata = translate(adata)
    assert mousified_adata.n_obs == adata.n_obs, "We lost cells during mapping, which should not happen!"
    assert mousified_adata.n_vars > 10000, "Very few genes (less than 10k) could be mapped! Expecting more!"

