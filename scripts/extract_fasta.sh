#!/bin/bash

# extract fasta sequences compressed in a tar.xz file
# path_tar=$1
# path_file=$2
outdir=data/reference_sequences

ls $outdir/lists_by_tar/ | while read f; \
 do tar -xvf data/batches_bacteria/$(basename $f .txt).tar.xz -C $outdir $(cat $outdir/lists_by_tar/$f | tr '\n' ' '); \
 done
