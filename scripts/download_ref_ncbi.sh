#!/usr/bin/bash

OUTDIR=data/reference_ncbi        # to download reference sequences in zip format
REFSEQS=data/ref-seqs             # to extract fasta (.fna) of all reference sequences
mkdir -p $OUTDIR
mkdir -p $REFSEQS
# ACCESSION=GCF_000006945.2

# download zip
tail -n +2 data/ncbi_reference_sequences.txt | cut -f2 | while read ACCESSION; do \
    datasets download genome accession $ACCESSION --include genome --filename $OUTDIR/$ACCESSION.zip; \
    sleep 0.01; \
        done

# unzip and get fasta files 
ls $OUTDIR | while read f; do \
    unzip $OUTDIR/$f -d $REFSEQS/$(basename $f .zip); \
        done

# get path to each reference sequence by assembly accession
ls $REFSEQS/*/*/*/*/*.fna | awk -F'/' '{print $3,$0}' > data/path_reference_by_accession.txt
