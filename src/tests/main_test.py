from scanpy import read  # will this fail if scanpy is not in requirements?
from mousipy import translate


def test_pancreas_biomart():
    # download mouse scRNA-seq dataset
    url_datadir = "https://github.com/theislab/scvelo_notebooks/raw/master/"
    url = f"{url_datadir}data/Pancreas/endocrinogenesis_day15.h5ad"
    adata = read("data/Pancreas/endocrinogenesis_day15.h5ad", backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()

    humanized_adata = translate(adata, source='biomart')
    assert humanized_adata.n_obs == adata.n_obs, "We lost cells during mapping, which should not happen!"
    assert humanized_adata.n_vars > 10000, "Very few genes (less than 10k) could be mapped! Expecting around 17.5k."


def test_human_PBMC_biomart():
    # download human scRNA-seq dataset
    url = 'http://falexwolf.de/data/pbmc3k_raw.h5ad'
    adata = read("data/Pancreas/pbmc3k_raw.h5ad", backup_url=url)
    adata.var_names_make_unique()

    mousified_adata = translate(adata, source='biomart')
    assert mousified_adata.n_obs == adata.n_obs, "We lost cells during mapping, which should not happen!"
    assert mousified_adata.n_vars > 10000, "Very few genes (less than 10k) could be mapped! Expecting more!"

def test_pancreas_hcop():
    # download mouse scRNA-seq dataset
    url_datadir = "https://github.com/theislab/scvelo_notebooks/raw/master/"
    url = f"{url_datadir}data/Pancreas/endocrinogenesis_day15.h5ad"
    adata = read("data/Pancreas/endocrinogenesis_day15.h5ad", backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()

    humanized_adata = translate(adata, source='hcop')
    assert humanized_adata.n_obs == adata.n_obs, "We lost cells during mapping, which should not happen!"
    assert humanized_adata.n_vars > 10000, "Very few genes (less than 10k) could be mapped! Expecting around 17.5k."


def test_PBMC_hcop():
    # download human scRNA-seq dataset
    url = 'http://falexwolf.de/data/pbmc3k_raw.h5ad'
    adata = read("data/Pancreas/pbmc3k_raw.h5ad", backup_url=url)
    adata.var_names_make_unique()

    mousified_adata = translate(adata, source='hcop')
    assert mousified_adata.n_obs == adata.n_obs, "We lost cells during mapping, which should not happen!"
    assert mousified_adata.n_vars > 10000, "Very few genes (less than 10k) could be mapped! Expecting more!"