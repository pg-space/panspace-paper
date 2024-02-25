.PHONY: all clean download_gold_standard download_bacteria fcgr create query config test_download_bacteria test_fcgr

SHELL=/usr/bin/env bash -eo pipefail
DATETIME=$(shell date -u +"%Y_%m_%dT%H_%M_%S")

.SECONDARY:

.SUFFIXES:

## config params
THREADS=$(shell grep "^threads:" pipeline/config/params.yaml | awk '{print $$2}')
NVIDIA_GPU=$(shell grep "^nvidia_gpu:" pipeline/config/params.yaml | awk '{print $$2}')

## Snakemake params
SMK_PARAMS=--cores ${THREADS} --rerun-incomplete --printshellcmds --keep-going --use-conda
SMK_PARAMS_GPU=--cores ${THREADS} --rerun-incomplete --printshellcmds --keep-going --use-conda --resources nvidia_gpu=$(NVIDIA_GPU)

## General commands ##
all: ## Run everything (the default rule)
	echo "Hello there" 

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
	snakemake -s pipeline/Snakefile --until download_bacteria $(SMK_PARAMS) 
	snakemake -s pipeline/Snakefile --until fcgr $(SMK_PARAMS)

create:
	echo "TODO"

query:
	echo "TODO"

####################
##     Tests      ##
## test_<command> ##
####################

test_download_bacteria:
	snakemake -s pipeline/Snakefile --until download_bacteria $(SMK_PARAMS) --config subset=test

test_fcgr:
	snakemake -s pipeline/Snakefile --until download_bacteria $(SMK_PARAMS) --config subset=test
	snakemake -s pipeline/Snakefile --until fcgr $(SMK_PARAMS) --config subset=test

test_create:
	echo "TODO"

test_query:
	echo "TODO"

###############
## Reporting ##
###############

config: ## Print configuration without comments
	@cat pipeline/config/params.yaml \
		| perl -pe 's/ *#.*//g' \
		| grep --color='auto' -E '.*\:'
	@#| grep -Ev ^$$

# report: ## Generate Snakemake report
# 	snakemake -s pipeline/Snakefile --report
