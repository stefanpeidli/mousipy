![GitHub Workflow Status](https://img.shields.io/github/workflow/status/stefanpeidli/mousipy/Python%20package)
![GitHub issues](https://img.shields.io/github/issues-raw/stefanpeidli/mousipy)
![PyPI - Downloads](https://img.shields.io/pypi/dm/mousipy?label=pip%20downloads)
![PyPI](https://img.shields.io/pypi/v/mousipy?label=PyPI%20version)

# mousipy
A python package that translates an AnnData single cell object from scanpy with mouse gene symbols into one with human gene symbols by mapping orthologs from biomart.


# Why?
Many people just uppercase a mouse gene symbol to get the human ortholog in scRNA-seq data. This works in most cases, but fails for some.
For example, there is no Cd8b gene in mice since the correct mouse ortholog to the human gene CD8B is Cd8b1. The gene CD8B is a defining marker for CD8+ T cells
which would get lost by just uppercasing gene symbols but is correctly retained by mapping gene symbols with mousipy. Another example is CD16 (human gene FCGR3A), which has mouse ortholog Fcgr4.

# Install
Just install via pip:

```pip install mousipy```

# Usage example
```
import scvelo as scv
from mousipy import translate
adata = scv.datasets.pancreas()  # mouse scRNA-seq dataset
humanized_adata = translate(adata)
```

# How it works
In `mousipy/biomart` are lists of mouse (GRCm39) and human (GRCh38.p13) orthologs exported from [biomart](https://www.ensembl.org/biomart/).
First, for all mouse gene symbols in adata.var_names we check if there is an ortholog in these lists. Then, for each mouse gene
- if there is exactly one human ortholog, the gene symbol is translated directly
- if there is an entry for that gene in the list explicitly mapping it to no ortholog, it will be discarded
- if there are multiple different human orthologs, the gene's expression counts are added to **all** its orthologs
- if the gene is not found in the list, we make it uppercase (and hope that that is the ortholog)

# What is an ortholog?
Two genes in different species are called orthologs if they share a common ancestry. At some point in the past these genes must have underwent a specification event.
