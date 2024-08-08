configfile: "pipeline/config/params.yaml"
import json
from pathlib import Path
import datetime
import copy
from os.path import join as pjoin

DATADIR = config["datadir"]
KMER = config["kmer_size"]
OUTDIR=Path(config["outdir"])
NAME_EXPERIMENT=config["name_experiment"]

PATH_TRAIN = Path(OUTDIR).joinpath(f"{KMER}mer/{NAME_EXPERIMENT}/cross-validation")
BIN_FASTANI="/home/avila/tools/fastANI"

list_tarfiles = [p.stem for p in  Path(PATH_TRAIN).joinpath("confident-learning/lists-by-tar").glob("*txt")]

def get_flags(wildcards):
    return [ pjoin(PATH_TRAIN, "confident-learning","seqs-issues",f"{tarfile}.flag") for tarfile in list_tarfiles]
    

rule all:
    input:
        expand(
            pjoin(PATH_TRAIN, "confident-learning","seqs-issues","{tarfile}.flag"),
            tarfile=list_tarfiles),
        pjoin(PATH_TRAIN, "confident-learning","ani.tsv"),        


rule ani_on_tarfile:
    input: 
        tarfile=pjoin(DATADIR, "batches_bacteria", "{tarfile}.tar.xz"),
        list_extract=pjoin(PATH_TRAIN, "confident-learning" ,"lists-by-tar", "{tarfile}.txt"),
        metadata_issues = Path(PATH_TRAIN).joinpath("confident-learning/metadata_issues.tsv"),
    output:
        flag=pjoin(PATH_TRAIN, "confident-learning","seqs-issues","{tarfile}.flag"),
    log:
        pjoin(DATADIR, "logs", "decompress_tarxz-{tarfile}.log"),
    params:
        dir_sequences=pjoin(PATH_TRAIN,"confident-learning","seqs-issues"),
        dir_list_ani=PATH_TRAIN.joinpath("confident-learning","lists-ANI"),
        fastani=BIN_FASTANI,
        dirsave=PATH_TRAIN.joinpath("confident-learning","ani-results")
    resources:
        limit_space=1,
    #     disk_mb=20_000_000
    shell:        
        """
        tar -xvf {input.tarfile} -C {params.dir_sequences} $(cat {input.list_extract} | tr '\n' ' ')

        mkdir -p {params.dirsave}
        mkdir -p {params.dir_sequences}

        cat {input.list_extract} | while read f; do \
            {params.fastani} -q {params.dir_sequences}/$f --rl {params.dir_list_ani}/$(basename $f .fa).txt -o {params.dirsave}/$(basename $f .fa).txt; \
                done
        
        rm -r {params.dir_sequences}/{wildcards.tarfile}
        touch {output.flag}
        """

rule consolidate_ani:
    input:
        get_flags
    output:
        pjoin(PATH_TRAIN, "confident-learning","ani.tsv"),
    params:
        path_reference_by_accession = Path(DATADIR).joinpath("path_reference_by_accession.txt"),
        path_metadata_references = Path(DATADIR).joinpath("ncbi_reference_sequences.txt"), 
        path_cv = PATH_TRAIN,       
    conda:
        "panspace-cli"
    shell:
        "panspace data-curation consolidate-ani {params.path_cv} {params.path_reference_by_accession} {params.path_metadata_references}"