# panspace-paper

This repository contains code to reproduce the results of the paper
[insert name](www.thelink.com)

Using [panspace](www.linktorepo.com), we indexed the [661k bacterial dataset](www.linktopaper.com)
and evaluate the results on three gold standard datasets: GEBA, FDA-ARGOS, and NCTC3000

## 1. Create a virtual environment with snakemake
```
mamba env create -n snakemake -f pipeline/workflow/envs/smk.yaml
mamba activate snakemake
```
you can use `conda` / `miniconda` / `micromamba` instead of `mamba`.

___

## 2. Commands

```{bash}
# Data
make download_gold_standard             # GEBA, FDA-ARGOS and NTC3000 datasets
make download_bacteria                  # 661k bacterial dataset
make download_index_bacteria            # download index for 661k bacterial dataset TODO: upload to zenodo

# Query Index
make query                              # query 661k bacterial index

# Create Index   download->fcgr->cross-validation->index
make fcgr                               # create FCGR (it will download the data if needed)
make create_ae                          # create a new index: train a new encoder, create faiss-index
make create_ml                          # create a new index: train a new encoder, create faiss-index

# Utilities
make clean                              # delete all info in data/ directory

# Test
make test_download_bacteria
make test_fcgr
make test_create_ae
make test_create_ml
```

___

### 2.1 Query the index
We provide 2 indexes, one based in `autoencoder`, and the other one based on `metric-learning`. 
Choose between one of the following in the main configuration file `pipeline/config/params.yaml`
```yaml
index_model: metric-learning # Options: autoencoder, metric-learning 
```
default is `metric-learning`


*EXTENSIONS ACCEPTED*: `.fa` , `.fna`, `.fasta`


___

### How it works?
We work in the space of $k$-mer distributions of the assemblies
($k$ is fixed by the user). The **goal** is to compress each assembly 
into an $n$-dimensional vector, called embedding. 
These embeddings are used to create an index for the initial dataset,
A query to this index correspond to an input assembly, and the query result
is a list with the $N$-closests embeddings (representing indexed assemblies). 
Distance is computed with the Euclidean Distance. 

The embedding representation is created such that when querying the index, 
assemblies of the same species are retrieved.

- The $k$-mer distribution of each assembly is represented by its FCGR. 
- FCGRs are used to train an encoder (current implementations allow the use of autoencoders, or metric learning architectures)
- the encoder maps an FCGR to an n-dimensional vector, where each assembly is represented by another vector.
- a query consists on an (draft) assembly, and it returns a ranked list of indexed assemblies with the smaller euclidean distance between their vector representations.