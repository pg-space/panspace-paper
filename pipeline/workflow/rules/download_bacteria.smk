# !/bin/bash

DIRSAVE="data/bacteria"
BASEPATH="https://zenodo.org/records/4602622/files"

from os.path import join as pjoin

def aggregate_downloaded_batches(wildcards):
    dirsave = checkpoints.download_batches.get(**wildcards).output[0]
    return list(Path(dirsave).rglob("*tar.xz"))

rule download_verification: 
    input:
        "data/batches_661k_bacteria.flag"
        # aggregate_batches

rule download_list:
    output:
        "data/batches_661k_bacteria.txt"    
    params: 
        basepath=BASEPATH
    shell:
        """
        mkdir -p $(dirname {output})
        wget {params.basepath}/_md5.txt -O {output}
        """
    
checkpoint download_batches:
    output:
        directory("data/bacteria")
    input: 
        "data/batches_661k_bacteria.txt"
    params:
        dirsave=DIRSAVE,
        link_zenodo=BASEPATH,
    shell:
        """
        mkdir -p {output[0]}
        cut {input} -d" " -f3 | head -n 2 | \
        while read f; 
            do wget {params.link_zenodo}/$f -O {params.dirsave}/$f;
        done
        """


rule aggregate_download_batches:
    input:
        aggregate_downloaded_batches
    output:
        "data/batches_661k_bacteria.flag"
    shell:
        "touch {output}"