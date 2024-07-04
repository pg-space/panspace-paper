"""
Starting from a dataset of FCGR (.npy files) in folder (PATH_FCGR), 
this pipeline performs a k-fold cross-validation training of the autoencoder. 

For each k-fold, train and test sets (list of paths) are created, 
- the train set is used to train an autoencoder
- the train set is used to create an Index
- the test set is used to identify outliers and mislabeled assemblies
- the test set is used to query the Index 
"""

# wildcards: kfold

configfile: "pipeline/config/params.yaml"
import json
from pathlib import Path
import datetime
import copy

KMER = config["kmer_size"]
OUTDIR=Path(config["outdir"])
DATADIR=Path(config["datadir"])
PATH_FCGR = Path(DATADIR).joinpath(f"fcgr/{KMER}mer")
NAME_EXPERIMENT=config["name_experiment"]
PATH_TRAIN=Path(OUTDIR).joinpath(f"{KMER}mer/{NAME_EXPERIMENT}/cross-validation")
ARCHITECTURE = config["architecture"]
LATENT_DIM = config["latent_dim"]
KFOLD = config["kfold"]
KFOLDS = [x+1 for x in range(KFOLD)]
LABELS = config["labels"]
LOSS = config["loss"]
HIDDEN_ACTIVATION = config["hidden_activation"]
OUTPUT_ACTIVATION = config["output_activation"]

# save params used to run these pipeline

_params = config
# _params["datetime"] =  datetime.datetime.now()
_params["kmer_size"] = KMER

path_save_params = Path(PATH_TRAIN).joinpath("params.yaml")
path_save_params.parent.mkdir(exist_ok=True, parents=True)
with open(path_save_params, "w") as fp: 
    json.dump(_params, fp, indent=1)


# TODO: add wildcards: loss  hidden_activation  output_activation

def get_outputs(wildcards):

    outputs = []
    for kfold in KFOLDS: 

        outputs.extend(
            Path(PATH_TRAIN).joinpath(f"{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/test/query_results.csv")
            for loss, hidden_activation, output_activation in zip(LOSS, HIDDEN_ACTIVATION, OUTPUT_ACTIVATION)
        ),    
    return outputs

rule all:
    input:
        get_outputs

rule kfold_split:
    output:
        expand( Path(PATH_TRAIN).joinpath("train_{kfold}-fold.txt") , kfold=KFOLDS),
        expand( Path(PATH_TRAIN).joinpath("test_{kfold}-fold.txt") , kfold=KFOLDS),
    input:
        list(PATH_FCGR.rglob("*.npy"))
    params: 
        datadir=PATH_FCGR, 
        outdir=PATH_TRAIN,
        kfold=KFOLD,
        labels=LABELS
    log:
        Path(PATH_TRAIN).joinpath("logs/kfold_split.log")
    conda: 
        "../envs/panspace.yaml"
    shell:
        """/usr/bin/time -v panspace trainer split-data-cross-validation \
        --datadir {params.datadir} \
        --outdir {params.outdir} \
        --kfold {params.kfold} \
        --labels {params.labels} 2> {log}
        """

rule train:
    output:
        PATH_TRAIN.joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/checkpoints").joinpath(f"weights-{ARCHITECTURE}.keras"),
    input:
        Path(PATH_TRAIN).joinpath("train_{kfold}-fold.txt"),
    log:
        Path(PATH_TRAIN).joinpath("logs/train_{loss}-{hidden_activation}-{output_activation}-{kfold}-fold.log")
    conda:
        "../envs/panspace.yaml"
    resources:
        nvidia_gpu=1
    params:
        outdir=lambda w: PATH_TRAIN.joinpath(f"{w.loss}-{w.hidden_activation}-{w.output_activation}-{w.kfold}-fold"),
        autoencoder=config["architecture"],
        latent_dim=config["latent_dim"],
        kmer=config["kmer_size"],
        epochs=config["epochs"],
        batch_size=config["batch_size"],
        optimizer=config["optimizer"],
        patiente_early_stopping=config["patiente_early_stopping"],
        patiente_learning_rate=config["patiente_learning_rate"],
        train_size=config["train_size"],
        seed=config["seed"],
        loss=lambda wildcards: wildcards.loss,
        hidden_activation=lambda wildcards: wildcards.hidden_activation,
        output_activation=lambda wildcards: wildcards.output_activation,
        preprocessing=config["preprocessing"]
    shell:
        """/usr/bin/time -v panspace trainer train-autoencoder \
        --training-list {input} \
        --outdir {params.outdir} \
        --autoencoder {params.autoencoder} \
        --latent-dim {params.latent_dim} \
        --kmer {params.kmer} \
        --preprocessing {params.preprocessing} \
        --epochs {params.epochs} \
        --batch-size {params.batch_size} \
        --batch-normalization \
        --optimizer {params.optimizer} \
        --patiente-early-stopping {params.patiente_early_stopping} \
        --patiente-learning-rate {params.patiente_learning_rate} \
        --hidden-activation {params.hidden_activation} \
        --output-activation {params.output_activation} \
        --train-size {params.train_size} \
        --loss {params.loss} 2> {log}
        """

rule extract_encoder:
    output:
        Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/models/encoder.keras"),
    input:
        PATH_TRAIN.joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/checkpoints").joinpath(f"weights-{ARCHITECTURE}.keras"),
    params:
        dir_save=lambda w: Path(PATH_TRAIN).joinpath(f"{w.loss}-{w.hidden_activation}-{w.output_activation}-{w.kfold}-fold/models")
    log:
        Path(PATH_TRAIN).joinpath("logs/extract_encoder_{loss}-{hidden_activation}-{output_activation}-{kfold}-fold.log")
    conda: 
        "../envs/panspace.yaml"
    shell:
        "/usr/bin/time -v panspace trainer split-autoencoder --path-checkpoint {input} --dirsave {params.dir_save} --encoder-only 2> {log}"

rule create_index:
    output:
        index=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/faiss-embeddings/panspace.index"),
        embeddings=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/faiss-embeddings/embeddings.npy"),
        id_embeddings=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/faiss-embeddings/labels.json"),
    input:
        encoder=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/models/encoder.keras"),
        files_to_index=Path(PATH_TRAIN).joinpath("train_{kfold}-fold.txt"),
    params:
        latent_dim=LATENT_DIM,
        kmer_size=KMER,
    resources:
        nvidia_gpu=1
    log:
        Path(PATH_TRAIN).joinpath("logs/create_index_{loss}-{hidden_activation}-{output_activation}-{kfold}-fold.log")
    conda: 
        "../envs/panspace.yaml"
    shell:
        """/usr/bin/time -v panspace index create \
        --files-to-index {input.files_to_index} \
        --col-labels 1 \
        --path-encoder {input.encoder} \
        --path-index {output.index}\
        --latent-dim {params.latent_dim} \
        --kmer-size {params.kmer_size} 2> {log}"""

rule test_index:
    output:
        embeddings=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/test/embeddings.npy"),
        query=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/test/query_results.csv"),
    input:
        path_index=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/faiss-embeddings/panspace.index"),
        path_encoder=Path(PATH_TRAIN).joinpath("{loss}-{hidden_activation}-{output_activation}-{kfold}-fold/models/encoder.keras"),
        files_to_query=Path(PATH_TRAIN).joinpath("test_{kfold}-fold.txt")
    params:
        outdir=lambda w: Path(PATH_TRAIN).joinpath(f"{w.loss}-{w.hidden_activation}-{w.output_activation}-{w.kfold}-fold/test"),
        kmer_size=KMER,
    resources:
        nvidia_gpu=1
    log:
        Path(PATH_TRAIN).joinpath("logs/test_index_{loss}-{hidden_activation}-{output_activation}-{kfold}-fold.log")
    conda: 
        "../envs/panspace.yaml"
    shell:        
        """
        /usr/bin/time -v panspace index query \
        --path-fcgr {input.files_to_query} \
        --path-encoder {input.path_encoder} \
        --path-index {input.path_index} \
        --col-labels 1 \
        --outdir {params.outdir} \
        --kmer-size {params.kmer_size} 2> {log}
        """