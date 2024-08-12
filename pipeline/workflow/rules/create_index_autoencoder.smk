"""
This script train a model with the modified data (remove outliers, fix labels issues)
then creates a faiss index with all the resulting data
"""

configfile: "pipeline/config/params.yaml"

from pathlib import Path

KMER = config["kmer_size"]
OUTDIR=config["outdir"]
DATADIR=Path(config["datadir"])
PATH_FCGR = Path(DATADIR).joinpath(f"fcgr/{KMER}mer")
NAME_EXPERIMENT=config["name_experiment"]
PATH_TRAIN=Path(OUTDIR).joinpath(f"{KMER}mer/{NAME_EXPERIMENT}")
ARCHITECTURE = config["architecture"]
LATENT_DIM = config["latent_dim"]
LABELS = config["labels"]
LOSS = config["loss"]
HIDDEN_ACTIVATION = config["hidden_activation"]
OUTPUT_ACTIVATION = config["output_activation"]


rule: 
    input: 
        Path(PATH_TRAIN).joinpath("faiss-embeddings/panspace.index")
        
# Generate a txt file with numpy files to use (first column) (remove outliers)
# and their labels (second colum) (modified labels included)
rule modify_labels_and_remove_outliers:
    input:
        # list(PATH_FCGR.rglob("*/*.npy")),
        outliers=PATH_TRAIN.joinpath("cross-validation/outliers/outliers.csv"),
        ani = PATH_TRAIN.joinpath("cross-validation/confident-learning/ani.tsv"), 
    output:
        PATH_TRAIN.joinpath("train_list.txt"),
    params:
        path_fcgr = PATH_FCGR,
        labels = LABELS,
    conda:
        "../envs/panspace.yaml"
    shell:
        "panspace data-curation remove-outliers-fix-labels {params.path_fcgr} {input.outliers} {input.ani} {params.labels} {output}"

rule train:
    output:
        PATH_TRAIN.joinpath(f"checkpoints/weights-{ARCHITECTURE}.keras"),
    input:
        training_list=PATH_TRAIN.joinpath("train_list.txt"),
    log:
        Path(PATH_TRAIN).joinpath("logs/train.log")
    conda:
        "../envs/panspace.yaml"
    resources:
        nvidia_gpu=1
    params:
        outdir=PATH_TRAIN,
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
    shell:
        """/usr/bin/time -v panspace trainer train-autoencoder \
        --training-list {input.training_list} \
        --outdir {params.outdir} \
        --autoencoder {params.autoencoder} \
        --latent-dim {params.latent_dim} \
        --kmer {params.kmer} \
        --epochs {params.epochs} \
        --batch-size {params.batch_size} \
        --optimizer {params.optimizer} \
        --patiente-early-stopping {params.patiente_early_stopping} \
        --patiente-learning-rate {params.patiente_learning_rate} \
        --train-size {params.train_size} \
        --seed {params.seed} 2> {log}
        """

rule encoder_decoder:
    output:
        Path(PATH_TRAIN).joinpath(f"models/encoder.keras"),
        Path(PATH_TRAIN).joinpath(f"models/decoder.keras")
    input:
        Path(PATH_TRAIN).joinpath(f"checkpoints/weights-{ARCHITECTURE}.keras")
    params:
        dir_save=Path(PATH_TRAIN).joinpath("models")
    log:
        Path(PATH_TRAIN).joinpath("logs/encoder_decoder.log")
    conda: 
        "../envs/panspace.yaml"
    shell:
        "/usr/bin/time -v panspace trainer split-autoencoder --path-checkpoint {input} --dirsave {params.dir_save} 2> {log}"

rule create_index:
    output:
        Path(PATH_TRAIN).joinpath("faiss-embeddings/panspace.index"),
        Path(PATH_TRAIN).joinpath("faiss-embeddings/embeddings.npy"),
        Path(PATH_TRAIN).joinpath("faiss-embeddings/id_embeddings.json"),
    input:
        Path(PATH_TRAIN).joinpath("models/encoder.keras"),
        PATH_TRAIN.joinpath("train_list.txt"),
    params:
        latent_dim=LATENT_DIM
    resources:
        nvidia_gpu=1
    log:
        Path(PATH_TRAIN).joinpath("logs/create_index.log")
    conda: 
        "../envs/panspace.yaml"
    shell:
        """/usr/bin/time -v panspace index create \
        --files-to-index {input[1]} \
        --col-labels 1 \
        --path-encoder {input[0]} \
        --path-index {output[0]}\
        --latent-dim {params.latent_dim} 2> {log}"""