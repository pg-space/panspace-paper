"""Download datasets
- GEBA
- FDA-ARGOS
- NCTC3000
"""

from pathlib import Path

PATH_DATASET_GOLD_STANDARD = Path("data/gold_standard")
PATH_DATASET_GOLD_STANDARD.mkdir(exist_ok=True, parents=True)

ACCESSIONS_NCBI={
    # "NCTC3000": "PRJEB6403",
    "GEBA": "PRJNA30815",
    # "FDA-ARGOS": "PRJNA231221",
}
datasets = list(ACCESSIONS_NCBI.keys())
print(datasets)

rule download_datasets:
    input: 
        # expand(
        #     PATH_DATASET_GOLD_STANDARD.joinpath("{dataset}/ncbi_dataset.zip"), dataset=datasets,    
        # ),
        PATH_DATASET_GOLD_STANDARD.joinpath("{dataset}/list_paths.txt"),

rule download:
    output:
        PATH_DATASET_GOLD_STANDARD.joinpath("{dataset}/ncbi_dataset.zip"),
    # input:
    #     pass
    params:
        accession=lambda w: ACCESSIONS_NCBI[w.dataset],
        filename=lambda w: PATH_DATASET_GOLD_STANDARD.joinpath(f"{w.dataset}/ncbi_dataset.zip"),
    conda: 
        "../envs/ncbi.yaml"
    shell:
        "/usr/bin/time -v datasets download genome accession {params.accession} --filename {params.filename}"


# # TODO: include rules below
checkpoint unzip:
    output:
        directory(PATH_DATASET_GOLD_STANDARD.joinpath("{dataset}/ncbi_dataset/data")) 
    input:
        PATH_DATASET_GOLD_STANDARD.joinpath("{dataset}/ncbi_dataset.zip"),
    params:
        outdir_zip=lambda w: PATH_DATASET_GOLD_STANDARD.joinpath(f"{w.dataset}"),
        path_datasets=PATH_DATASET_GOLD_STANDARD
    log:
        "logs/{dataset}-unzip.log"
    shell:
        """
        unzip {input} -d {params.outdir_zip}
        """        

def get_fasta(wildcards):
    path_dataset = Path(checkpoints.unzip.get(**wildcards).output[0])
    list_paths = path_dataset.rglob("*fna")
    print(list_paths)
    return list_paths


rule create_list_paths:
    output:
        PATH_DATASET_GOLD_STANDARD.joinpath("{dataset}/list_paths.txt"),
    input:
        get_fasta 
    params:
        outdir_zip=lambda w: PATH_DATASET_GOLD_STANDARD.joinpath(f"{w.dataset}"),
        path_datasets=PATH_DATASET_GOLD_STANDARD
    log:
        "logs/{dataset}-create_list_paths.log"
    shell:
        """
        ls {params.outdir_zip}/ncbi_dataset/data/*/*fna | while read f; do echo $f >> {output}; done
        """
