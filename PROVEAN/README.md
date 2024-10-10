# AcrIIa4 mutant predictions with PROVEAN

Specifically requires
- `ncbi-blast-2.4.0+`
- `provean-1.1.5` extracted from [`Sourceforge`](https://sourceforge.net/projects/provean/) download

```bash
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ncbi-blast-2.4.0+-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.4.0+-x64-linux.tar.gz
```
module load cd-hit
module load blast


## Install `PROVEAN` in local directory (no `sudo` required)
```bash
tar -zxvf provean-1.1.5.tar.gz
cd provean-1.1.5

# Install 
./configure --prefix=/data/${USER}/provean

make
make install
```

## Download v4 `nr` database from NCBI

```bash
# Retrieve v4 NR database
mkdir db 
cd db
wget --recursive \
    -e robots=off \
    --reject "indexl.html" \
    --no-host-directories \
    --cut-dirs=6 https://ftp.ncbi.nlm.nih.gov/blast/db/v4/ \
    -A "nr_v4*"

# Extract files one at a time; they share `nr.pal` and will collide if done simultaneously
for file in *.tar.gz; do
    tar -zxvf ${file}
done
```

## Edit PROVEAN configuration
```bash
vi /data/${USER}/provean/bin/provean.sh
```
| Variable | Defined as |
| ---------|------------|
| `BLAST_DB` | the `nr` database downloaded with `wget` from NCBI |
| `PSIBLAST` | the version within the extracted `ncbi-blast-2.4.0+` directory |
| `BLASTDBCMD` | the version within the extracted `ncbi-blast-2.4.0+` directory |
| `CD_HIT` | the version pre-included on Biowulf, `/usr/local/apps/cd-hit/cdhit-4.8.1/cd-hit` |

## Build list of HGVS mutants

[`make-hgvs-file.py`](make-hgvs-file.py) generates [`hgvs-muts.txt`](hgvs-muts.txt)


## Run PROVEAN

```bash
bash /data/${USER}/provean/bin/provean.sh \
    --query AcrIIa4.fa \
    --variation hgvs-muts.txt \
    --save_supporting_set AcrIIa4.sss > provean-out.txt
```
