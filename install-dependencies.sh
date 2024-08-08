#!/bin/bash

set -e

# Install the conda environment:
mamba env create -f requirements.yaml --force

source $(conda info --base)/etc/profile.d/conda.sh
# Since the name changes with each version, make sure to pull the version
# name out of the yml file
ENV_NAME=$(awk '/name: mobilome_submit_v/{print $NF}' requirements.yaml)
echo $ENV_NAME
conda activate $ENV_NAME

# Install non-conda packages (we use the cloud-0 mirror so that the download
# will be fast in any part of the world)
#Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")'
#Rscript -e 'BiocManager::install("ggtree")'
Rscript -e 'devtools::install_github("thackl/thacklr")'
Rscript -e 'devtools::install_github("thackl/gggenomes")'

# Deactivate the conda environment
conda deactivate
