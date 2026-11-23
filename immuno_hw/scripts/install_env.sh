#!/bin/bash
#SBATCH --job-name=install_mixcr
#SBATCH --partition=AMD9554-common
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=outputs/install_env_%j.log

source /home/STUDY/FBMF/studfbmf02_13/miniforge3/etc/profile.d/conda.sh

mamba env create -f environment.yml -y
conda activate immune-mixcr
python -m ipykernel install --user --name immune-mixcr --display-name "Python (immune-mixcr)"
mixcr --version
