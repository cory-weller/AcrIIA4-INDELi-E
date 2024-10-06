#!/usr/bin/env bash
#SBATCH --mem 10G
#SBATCH --time 0:10:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --partition quick

cd /data/SBGE/cory/AcrIIA4-vep/esm-variants

python3 esm_score_missense_mutations.py \
    --input-fasta-file ../AcrIIa4.fa \
     --output-csv-file ../esm-missense.csv