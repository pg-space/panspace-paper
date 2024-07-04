configfile: "pipeline/config/params.yaml"
import json
from pathlib import Path
import datetime
import copy

KMER = config["kmer_size"]
OUTDIR=Path(config["outdir"])
NAME_EXPERIMENT=config["name_experiment"]

PATH_TRAIN = Path(OUTDIR).joinpath(f"{KMER}mer/{NAME_EXPERIMENT}/cross-validation")
PERCENTIL = config["outliers"]["percentile_avg_distance"] 

KFOLD = config["kfold"]
KFOLDS = [x+1 for x in range(KFOLD)]

LOSS = config["loss"]
HIDDEN_ACTIVATION = config["hidden_activation"]
OUTPUT_ACTIVATION = config["output_activation"]


def get_outputs(wildcards):

    outputs = []
    for kfold in KFOLDS: 

        outputs.extend(
            Path(PATH_TRAIN).joinpath(f"{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/ood/outliers_avg_dist_percentile.csv")
            for loss, hidden_activation, output_activation in zip(LOSS, HIDDEN_ACTIVATION, OUTPUT_ACTIVATION)
        )
        
    return outputs

rule all:
    input:
        get_outputs,
        Path(PATH_TRAIN).joinpath(f"outliers/outliers.csv")



rule outlier_detection:
    output:
        path_outliers = Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/ood/outliers_avg_dist_percentile.csv"),
    input:
        path_index=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/faiss-embeddings/panspace.index"),
        train_embeddings=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/faiss-embeddings/embeddings.npy"),
        test_embeddings=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/test/embeddings.npy"),
    params:
        # train_metadata=lambda w: Path(PATH_TRAIN).joinpath(f"train_{w.kfold}-fold.txt"),
        test_metadata=lambda w: Path(PATH_TRAIN).joinpath(f"test_{w.kfold}-fold.txt"),
        outdir=lambda w: Path(PATH_TRAIN).joinpath(f"{w.loss}-{w.hidden_activation}-{w.output_activation}-{w.kfold}-fold/ood"),
        neighbors=10,
        percentile=PERCENTIL,
    conda: 
        "../envs/panspace.yaml"
    log:
        log=Path(PATH_TRAIN).joinpath("logs/outlier_detection_{loss}-{hidden_activation}-{output_activation}-{kfold}-fold.log")
    shell:
        """
        /usr/bin/time -v panspace data-curation ood \
        --path-index {input.path_index} \
        --path-train-embeddings {input.train_embeddings} \
        --path-test-embeddings {input.test_embeddings}\
        --path-test-metadata {params.test_metadata} \
        --outdir {params.outdir} \
        --neighbors {params.neighbors} \
        --percentile {params.percentile} 2> {log} 
        """

rule merge_outliers:
    output:
        path_outliers = Path(PATH_TRAIN).joinpath(f"outliers/outliers.csv"),
    input:
        get_outputs
    conda: 
        "../envs/panspace.yaml"
    log:
        log=Path(PATH_TRAIN).joinpath("logs/merge_outliers.log")
    shell:
        """
        /usr/bin/time -v panspace data-curation utils-join-df \
        {input} -o {output} 2> {log}
        """