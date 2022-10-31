import scvelo as scv
from mousipy import translate

def test_pancreas():
    adata = scv.datasets.pancreas()  # mouse scRNA-seq dataset
    humanized_adata = translate(adata)
    assert humanized_adata.n_obs == adata.n_obs, "We lost cells during mapping, which should not happen!"
    assert humanized_adata.n_vars > 10000, "Very few genes (less than 10k) could be mapped! Expecting around 17.5k."
