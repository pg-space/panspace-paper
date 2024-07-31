#!/usr/bin/bash

BIN_FASTANI=/home/avila/tools/fastANI 
dirtarfiles=data/batches_bacteria
cfdir=experiments-paper/6mer/07_25_2024-autoencoder/cross-validation/confident-learning
outdir=$cfdir/seqs-issues
dir_list_ani=$cfdir/lists-ANI
dir_list_by_tar=$cfdir/lists-by-tar
tarfile=burkholderia_pseudomallei__01 #acinetobacter_baumannii__01
# results=
dirsave=$outdir/ani-results
mkdir -p $dirsave

tar -xvf $dirtarfiles/$tarfile.tar.xz -C $outdir $(cat $dir_list_by_tar/$tarfile.txt | tr '\n' ' ')

cat $dir_list_by_tar/$tarfile.txt | while read f; do \
    $BIN_FASTANI -q $outdir/$f --rl $dir_list_ani/$(basename $f .fa).txt -o $dirsave/$(basename $f .fa).txt --minFraction 0.001; \
        done
# rm -r $outdir/$tarfile