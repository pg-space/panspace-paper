.PHONY: 
	download_gold_standard
	download_bacteria
	fcgr
	create
	query
	test_download_bacteria
	test_fcgr
	test_create
	test_query
	config
	clean

SHELL=/usr/bin/env bash -eo pipefail
DATETIME=$(shell date -u +"%Y_%m_%dT%H_%M_%S")

.SECONDARY:

.SUFFIXES:

## config params
THREADS=$(shell grep "^threads:" pipeline/config/params.yaml | awk '{print $$2}')
NVIDIA_GPU=$(shell grep "^nvidia_gpu:" pipeline/config/params.yaml | awk '{print $$2}')

## Snakemake params
SMK_PARAMS=--cores ${THREADS} --rerun-incomplete --printshellcmds --keep-going --use-conda
SMK_PARAMS_GPU=--cores ${THREADS} --rerun-incomplete --printshellcmds --keep-going --use-conda --resources nvidia_gpu=${NVIDIA_GPU}

## --- Commands --- 

####################
## Download data ##
####################

download_gold_standard: ## Download assemblies from NCBI: GEBA, FDA-ARGOS, NCTC3000
	snakemake -s pipeline/Snakefile --until download_gold_standard $(SMK_PARAMS)

download_bacteria: ## Download 661k bacterial assembly dataset
	snakemake -s pipeline/Snakefile --until download_bacteria $(SMK_PARAMS) 

####################
## Pipeline steps ##
####################

fcgr:
	make download_bacteria
	snakemake -s pipeline/Snakefile --until fcgr $(SMK_PARAMS)

# metric learning
create_ml:
	make fcgr 
	snakemake -s pipeline/workflow/rules/cross_validation_metric_learning.smk $(SMK_PARAMS_GPU)

# autoencoder
create_ae:
	make fcgr 
	snakemake -s pipeline/workflow/rules/cross_validation_autoencoder.smk $(SMK_PARAMS_GPU)

query:
	echo "TODO"


####################
##     Utils      ##
####################

config: ## Print configuration without comments
	@cat pipeline/config/params.yaml \
		| perl -pe 's/ *#.*//g' \
		| grep --color='auto' -E '.*\:'
	@#| grep -Ev ^$$

clean:
	rm -rf data/bacteria_test
	rm -rf data/bacteria_all
	rm -rf data/gold_standard
	rm -rf data/kmer-count 
	rm -rf data/fcgr
	rm -rf data/assembly
	rm -rf data/logs
	rm -rf data/list*txt
	rm -rf data/*flag

clean_output:
	rm -rf output/*
	
####################
##     Tests      ##
## test_<command> ##
####################

test_download_bacteria:
	snakemake -s pipeline/Snakefile --until download_bacteria $(SMK_PARAMS) --config subset=test

test_fcgr:
	make test_download_bacteria
	snakemake -s pipeline/Snakefile --until fcgr $(SMK_PARAMS) --config subset=test

test_create_ml:
	make test_fcgr
	snakemake -s pipeline/workflow/rules/cross_validation_metric_learning.smk $(SMK_PARAMS_GPU) --config subset=test architecture=CNNFCGR outdir=output/test-metric-learning hidden_activation=[relu] loss=[triplet_semihard_loss] 

test_create_ae:
	make test_fcgr 
	snakemake -s pipeline/workflow/rules/cross_validation_autoencoder.smk $(SMK_PARAMS_GPU) --config subset=test architecture=AutoencoderFCGR outdir=output/test-autoencoder hidden_activation=[relu] output_activation=[sigmoid] loss=[binary_crossentropy]
	
test_query:
	echo "TODO"

