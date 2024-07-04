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

print(config)

def get_outputs(wildcards):

    outputs = []
    for kfold in KFOLDS: 

        outputs.extend(
            Path(PATH_TRAIN).joinpath(f"{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/confident-learning/pred_probs.npy")
            for loss, hidden_activation, output_activation in zip(LOSS, HIDDEN_ACTIVATION, OUTPUT_ACTIVATION)
        )
        
    return outputs


def get_outputs_labels(wildcards):

    outputs = []
    for kfold in KFOLDS: 

        outputs.extend(
            Path(PATH_TRAIN).joinpath(f"{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/confident-learning/labels.npy")
            for loss, hidden_activation, output_activation in zip(LOSS, HIDDEN_ACTIVATION, OUTPUT_ACTIVATION)
        )
        
    return outputs

rule all:
    input:
        # get_outputs,
        path_pred_probs = Path(PATH_TRAIN).joinpath(f"confident-learning/pred_probs.npy"),
        path_labels = Path(PATH_TRAIN).joinpath(f"confident-learning/labels.npy"),

rule pred_probs:
    output:
        path_pred_probs = Path(PATH_TRAIN).joinpath(f"{{loss}}-{{hidden_activation}}-{{output_activation}}-{{kfold}}-fold/confident-learning/pred_probs.npy"),
        path_label_gt = Path(PATH_TRAIN).joinpath(f"{{loss}}-{{hidden_activation}}-{{output_activation}}-{{kfold}}-fold/confident-learning/labels.npy"),
        # path_pred_labels = Path(PATH_TRAIN).joinpath(f"{{loss}}-{{hidden_activation}}-{{output_activation}}-{{kfold}}-fold/confident_learning/pred_labels.csv"),
    input:
        train_embeddings=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/faiss-embeddings/embeddings.npy"),
        test_embeddings=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/test/embeddings.npy"),
        # train_labels=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/faiss-embeddings/labels.json"),
        # test_labels=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/test/labels.json"),
    params:
        train_labels=lambda w: Path(PATH_TRAIN).joinpath(f"train_{w.kfold}-fold.txt"),
        test_labels=lambda w: Path(PATH_TRAIN).joinpath(f"test_{w.kfold}-fold.txt"),
        outdir=lambda w: Path(PATH_TRAIN).joinpath(f"{w.loss}-{w.hidden_activation}-{w.output_activation}-{w.kfold}-fold/confident-learning"),
    conda: 
        "../envs/panspace.yaml"
    log:
        log=Path(PATH_TRAIN).joinpath("logs/confident_learning_{loss}-{hidden_activation}-{output_activation}-{kfold}-fold.log")
    shell:
        """
        /usr/bin/time -v panspace data-curation pred-scores \
        --path-train-embeddings {input.train_embeddings} \
        --path-test-embeddings {input.test_embeddings} \
        --path-train-labels {params.train_labels} \
        --path-test-labels {params.test_labels} \
        --outdir {params.outdir} 2> {log} 
        """


rule merge_pred_probs:
    output:
        path_pred_probs = Path(PATH_TRAIN).joinpath(f"confident-learning/pred_probs.npy"),
    input:
        get_outputs
    params:
        path_save=lambda w: Path(PATH_TRAIN).joinpath("confident-learning/pred_probs.npy"),
    conda: 
        "../envs/panspace.yaml"
    log:
        log=Path(PATH_TRAIN).joinpath("logs/confident_learning.log")
    shell:
        """
        /usr/bin/time -v panspace data-curation utils-join-npy \
        {input} -o {params.path_save} 2> {log}
        """


rule merge_labels:
    output:
        path_labels = Path(PATH_TRAIN).joinpath(f"confident-learning/labels.npy"),
    input:
        get_outputs_labels
    params:
        path_save=lambda w: Path(PATH_TRAIN).joinpath("confident-learning/labels.npy"),
    conda: 
        "../envs/panspace.yaml"
    log:
        log=Path(PATH_TRAIN).joinpath("logs/confident_learning.log")
    shell:
        """
        /usr/bin/time -v panspace data-curation utils-join-npy \
        {input} -o {params.path_save} 2> {log}
        """