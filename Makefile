.PHONY: all conda test fcgr download_gold_standard download_bacteria fcgr_bacteria fcgr_gold_standard config report

SHELL=/usr/bin/env bash -eo pipefail
DATETIME=$(shell date -u +"%Y_%m_%dT%H_%M_%S")

.SECONDARY:

.SUFFIXES:

THREADS=$(shell grep "^threads:" pipeline/config/params.yaml | awk '{print $$2}')
# MAX_DOWNLOAD_THREADS=$(shell grep "^max_download_threads"  pipeline/config/params.yaml | awk '{print $$2}')
# # DOWNLOAD_RETRIES=$(shell grep "^download_retries" config.yaml | awk '{print $$2}')
# MAX_IO_HEAVY_THREADS=$(shell grep "^max_io_heavy_threads" pipeline/config/params.yaml | awk '{print $$2}')
# MAX_RAM_MB=$(shell grep "^max_ram_gb:" pipeline/config/params.yaml | awk '{print $$2*1024}')

# resources
NVIDIA_GPU=$(shell grep "^nvidia_gpu:" pipeline/config/params.yaml | awk '{print $$2}')

SMK_PARAMS=--cores ${THREADS} --rerun-incomplete --printshellcmds --keep-going --use-conda
SMK_PARAMS_GPU=--cores ${THREADS} --rerun-incomplete --printshellcmds --keep-going --use-conda --resources nvidia_gpu=$(NVIDIA_GPU)

# ifeq ($(SMK_CLUSTER_ARGS),)
#     # configure local runSHELL
#     SMK_PARAMS=--cores ${THREADS} --rerun-incomplete --printshellcmds --keep-going --use-conda --resources max_download_threads=$(MAX_DOWNLOAD_THREADS) max_io_heavy_threads=$(MAX_IO_HEAVY_THREADS) max_ram_mb=$(MAX_RAM_MB)
# else
#     # configure cluster run
#     SMK_PARAMS=--cores all --rerun-incomplete --printshellcmds --keep-going --use-conda --resources max_download_threads=10000000 max_io_heavy_threads=10000000 max_ram_mb=1000000000 $(SMK_CLUSTER_ARGS)
# endif

# DOWNLOAD_PARAMS=--cores $(MAX_DOWNLOAD_THREADS) -j $(MAX_DOWNLOAD_THREADS) --restart-times $(DOWNLOAD_RETRIES)


######################
## General commands ##
######################
all: ## Run everything (the default rule)
	echo "Hello there" 
	# make download_gold_standard
	# make download_bacteria
	# make create_index
	# make query

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
####################

test_download_bacteria:
	snakemake -s pipeline/Snakefile --until download_bacteria $(SMK_PARAMS) --config subset=test

test_fcgr:
	snakemake -s pipeline/Snakefile --until download_bacteria $(SMK_PARAMS) --config subset=test
	snakemake -s pipeline/Snakefile --until fcgr $(SMK_PARAMS) --config subset=test

# test:


# conda: ## Create the conda environments
# 	snakemake -s pipeline/Snakefile $(SMK_PARAMS) --conda-create-envs-only

# download_gold_standard: ## Download assemblies from NCBI: GEBA, FDA-ARGOS, NCTC3000
# 	snakemake -s pipeline/Snakefile --until download_gold_standard $(SMK_PARAMS) $(DOWNLOAD_PARAMS)

# download: ## Download only the assemblies
# 	snakemake -s pipeline/Snakefile --until download_bacteria $(SMK_PARAMS) $(DOWNLOAD_PARAMS)

# download_cobs: ## Download only the COBS indexes
# 	snakemake download_cobs_batches $(SMK_PARAMS) $(DOWNLOAD_PARAMS)

# match: ## Match queries using COBS (queries -> candidates)
# 	scripts/benchmark.py --log logs/benchmarks/match_$(DATETIME).txt "snakemake match $(SMK_PARAMS)"

###############
## Reporting ##
###############

config: ## Print configuration without comments
	@cat config.yaml \
		| perl -pe 's/ *#.*//g' \
		| grep --color='auto' -E '.*\:'
	@#| grep -Ev ^$$

report: ## Generate Snakemake report
	snakemake --report



# ##########
# ## Misc ##
# ##########
# cluster_slurm: ## Submit to a SLURM cluster
# 	sbatch \
#         -c 10 \
#         --mem=80GB \
#         -t 0-08:00:00 \
#         --wrap="make"

# cluster_lsf_test: ## Submit the test pipeline to LSF cluster
# 	scripts/check_if_config_is_ok_for_cluster_run.py
# 	scripts/submit_lsf.sh test

# cluster_lsf: ## Submit to LSF cluster
# 	scripts/check_if_config_is_ok_for_cluster_run.py
# 	scripts/submit_lsf.sh

# format: ## Reformat Python and Snakemake files
# 	yapf -i */*.py
# 	snakefmt Snakefile
