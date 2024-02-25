configfile: "pipeline/config/params.yaml"
import json
import tarfile
from os.path import join as pjoin
from tqdm import tqdm
from collections import defaultdict
from pathlib import Path

# params
KMER_SIZE=config["kmer_size"]
OUTDIR=Path(config["outdir"])
SUBSET=config["subset"]

# # --- check all tarfiles ---
DIR_TARFILES=OUTDIR.joinpath(f"bacteria_{SUBSET}")
# TARFILES,= glob_wildcards(pjoin(DIR_TARFILES,"{tarfile}"+".tar.xz"))
# print(TARFILES)
# print(DIR_TARFILES)
# with open(OUTDIR.joinpath(f"batches_661k_bacteria_{SUBSET}.txt")) as fp:
### ---- FCGR ----
def load_batches(subset):
    list_files=[]
    with open(f"data/batches_661k_bacteria_{subset}.txt") as fp:
        
        for line in fp.readlines():
            tarxz = line.replace("\n","").strip().split(" ")[-1]
            name = tarxz.replace(".tar.xz","")
            list_files.append(name)
    return list_files

TARFILES=load_batches(SUBSET)

rule fcgr_verification:
    input:
        expand( pjoin(OUTDIR, "fcgr","{tarfile}_fcgr.flag"), tarfile=TARFILES)

# outut fasta files in assembly/ directory
checkpoint decompress_tarxz:
    input: 
        tarfile=pjoin(OUTDIR,f"bacteria_{SUBSET}", "{tarfile}" + ".tar.xz"),
        flag_downloaded=pjoin(OUTDIR,f"batches_661k_bacteria_{SUBSET}.flag"),
    output:
        directory(pjoin(OUTDIR, "assembly" ,"{tarfile}")),
    log:
        pjoin(OUTDIR, "logs", "decompress_tarxz-{tarfile}.log"),
    params:
        outdir=pjoin(OUTDIR,"assembly"),
    resources:
        limit_space=5,
    #     disk_mb=20_000_000
    shell:
        """
        mkdir -p {params.outdir}
        /usr/bin/time -v tar -xvf {input.tarfile} -C {params.outdir} 2> {log}
        """

# store KMC output (.kmc_pre and .kmc_suf) in fcgr/ directory
rule count_kmers:
    input:
        pjoin(OUTDIR, "assembly", "{tarfile}", "{fasta}.fa")
    output:
        pjoin(OUTDIR, "kmer-count","{tarfile}","{fasta}.kmc_pre"),
        pjoin(OUTDIR, "kmer-count","{tarfile}","{fasta}.kmc_suf"),
    log:
        pjoin(OUTDIR, "logs", "count_kmers-{tarfile}-{fasta}.log")
    params:
        kmer=KMER_SIZE,
        out=lambda w: pjoin(OUTDIR, "kmer-count",f"{w.tarfile}",f"{w.fasta}"),
    conda:
        "../envs/kmc.yaml"
    # resources:
        # limit_space=1,
        # disk_mb=20_000_000,
    # priority:
    #     100
    shell:
        """
        /usr/bin/time -v kmc -v -k{params.kmer} -m4 -sm -ci0 -cs65535 -b -t4 -fm {input} {params.out} . 2> {log}
        rm -r {input}
        """


def aggregate_fasta_kmc(wildcards,):
    "Helper function to collect all .kmc_suf files resulting from running KMC on the set of assemblies of a tarfile"
    
    output_tarfile = checkpoints.decompress_tarxz.get(**wildcards).output[0]
    list_fasta = glob_wildcards( pjoin(output_tarfile, "{fasta}.fa") ).fasta
    outdir = pjoin(OUTDIR, "kmer-count",f"{wildcards.tarfile}")
    return expand( pjoin(outdir,"{fasta}.kmc_suf"), fasta=list_fasta)    

rule list_path_fasta:
    input:  
        aggregate_fasta_kmc
    output: 
        pjoin(OUTDIR, "list_path_kmc_{tarfile}.txt")
    params:
        kmerdir=lambda w: pjoin(OUTDIR,"kmer-count",f"{w.tarfile}"),
        fcgrdir=lambda w: pjoin(OUTDIR,"fcgr",f"{w.tarfile}"),
        parent_fcgrdir = lambda w: pjoin(OUTDIR,"fcgr")
    log:
        OUTDIR.joinpath("logs/list_path_fasta-{tarfile}.log")
    shell:
        """
        /usr/bin/time -v ls {params.kmerdir}/*.kmc_suf | while read f; do echo ${{f::-8}} >> {output} ; done 2> {log}
        """


checkpoint save_fcgr_as_numpy:
    input:
        pjoin(OUTDIR, "list_path_kmc_{tarfile}.txt")
    output:
        directory(pjoin(OUTDIR,"fcgr","{tarfile}"))
    params:
        kmer=KMER_SIZE,
        kmerdir=lambda w: pjoin(OUTDIR,"kmer-count",f"{w.tarfile}"),
        fcgrdir=lambda w: pjoin(OUTDIR,"fcgr",f"{w.tarfile}"),
        bin_fcgr=config["bin_fcgr"]
    log:
        OUTDIR.joinpath("logs/fcgr-{tarfile}.log")
    priority:
        100
    conda: 
        "../envs/panspace.yaml"
    shell:
        """
        /usr/bin/time -v {params.bin_fcgr} {input} 2> {log}
        mkdir -p {params.fcgrdir} 
        mv {params.kmerdir}/*.npy {params.fcgrdir}
        """


def aggregate_numpy_fcgr(wildcards,):
    "Helper function to collect all FCGR .npy files generated for a set of assemblies of a tarfile"
    
    output_tarfile = checkpoints.save_fcgr_as_numpy.get(**wildcards).output[0]
    list_fasta = glob_wildcards( pjoin(output_tarfile, "{fasta}.fa") ).fasta
    return expand( pjoin(OUTDIR, "fcgr",f"{wildcards.tarfile}","{fasta}.npy"), fasta=list_fasta)    


rule fcgr_aggregate:
    input: 
        aggregate_numpy_fcgr
    output: 
        touch( pjoin(OUTDIR, "fcgr","{tarfile}_fcgr.flag"))
    priority:
        200
    params:     
        kmerdir=lambda w: pjoin(OUTDIR,"kmer-count",f"{w.tarfile}"),
    shell:
        "rm -rf {params.kmerdir}"
