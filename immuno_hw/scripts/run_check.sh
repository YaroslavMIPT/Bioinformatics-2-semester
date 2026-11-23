#!/bin/bash
#SBATCH --job-name=check_stats
#SBATCH --partition=AMD9554-common
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=outputs/check_stats_%j.log

source /home/STUDY/FBMF/studfbmf02_13/miniforge3/etc/profile.d/conda.sh
conda activate immune-mixcr

python scripts/check_stats.py
