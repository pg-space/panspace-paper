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

## 2. Query the index
or create a new one.

```{bash}
# Data
make download_gold_standard             # GEBA, FDA-ARGOS and NTC3000 datasets
make download_bacteria                  # 661k bacterial dataset
make download_index_bacteria            # download index for 661k bacterial dataset TODO: upload to zenodo

# Query Index
make query                              # query 661k bacterial index

# Create Index
make fcgr                               # create FCGR for 661k bacterial dataset
make create                             # create a new index: train a new encoder, create faiss-index 

# Utilities
make clean                              # delete all info in data/ directory
```

___


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
- The goal is to m