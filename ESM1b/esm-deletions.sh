#!/usr/bin/env bash
#SBATCH --mem 10G
#SBATCH --time 3:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --partition gpu
#SBATCH --gres gpu:v100x:1

module load CUDA

cd /data/SBGE/cory/AcrIIA4-vep/esm-variants

python3 esm_score_multi_residue_mutations.py \
    --input-csv-file ../muts.csv \
    --output-csv-file ../esm-deletions.csv