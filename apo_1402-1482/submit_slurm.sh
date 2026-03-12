#!/bin/bash
#SBATCH --job-name=apo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=05-00:00:00
#SBATCH --mem=16GB
#SBATCH --gpus=3
#SBATCH -o output.txt
#SBATCH -e error.txt
#SBATCH -p gpu

conda activate openmm
bash run.sh 2
