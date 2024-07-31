#!/usr/bin/bash

OUTDIR=data/reference_ncbi
mkdir -p $OUTDIR

# ACCESSION=GCF_000006945.2

tail -n +2 data/ncbi_reference_sequences.txt | cut -f2 | while read ACCESSION; do \
    datasets download genome accession $ACCESSION --include genome --filename $OUTDIR/$ACCESSION.zip; \
    sleep 0.01; \
        done