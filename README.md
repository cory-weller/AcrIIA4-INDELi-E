# AcrIIA4 mutant predictions

Prep peptide fasta
```bash
echo '>AcrIIa4' > AcrIIa4.fa
echo 'MNINDLIREIKNKDYTVKLSGTDSNSITQLIIRVNNDGNEYVISESENES
IVEKFISAFKNGWNQEYEDEEEFYNDMQTITLKSELN' >> AcrIIa4.fa
```

# Install requirements for variant prediction
```bash
pip3 install --user tqdm numpy pandas biopython torch fair-esm
git clone https://github.com/ntranoslab/esm-variants.git
```

# Prepare multiple-site deletion file
Requires a three-column file with `wt_seq`, `mut_seq`, `start_pos`.
Input built with [`build-muts.py`](build-muts.py)
```bash
python3 build-muts.py > muts.csv
```

# Run the models
```bash
sbatch esm-deletions.sh     # Needs GPU
bash esm-missense.sh   # Fast enough without GPU
```

# Plot results

Plots generated with [`plot-variants.R`](plot-variants.R)

![](esm-missense.png)
![](esm-deletions.png)