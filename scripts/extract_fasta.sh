#!/bin/bash

# extract fasta sequences compressed in a tar.xz file

outdir=experiments-paper/6mer/07_25_2024-autoencoder/cross-validation/confident-learning/seqs-issues
dir_list_by_tar=experiments-paper/6mer/07_25_2024-autoencoder/cross-validation/confident-learning/lists-by-tar

mkdir -p $outdir

ls $dir_list_by_tar | while read f; \
 do tar -xvf data/batches_bacteria/$(basename $f .txt).tar.xz -C $outdir $(cat $dir_list_by_tar/$f | tr '\n' ' '); \
 done
