

nrdb='/fdb/blastdb/nr'

wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ncbi-blast-2.4.0+-x64-linux.tar.gz

gunzip ncbi-blast-2.4.0+-x64-linux.tar.gz

module load cd-hit
module load blast

```bash
tar -zxvf provean-1.1.5.tar.gz
cd provean-1.1.5


./configure BLAST_DB=/fdb/blastdb/nr
./configure --prefix=/data/wellerca/provean

module load blast
module load cd-hit

./configure PSIBLAST=/usr/local/apps/blast/ncbi-blast-2.15.0+/bin/psiblast
./configure CDHIT=/usr/local/apps/cd-hit/cdhit-4.8.1/cd-hit

make
make install

bash /data/wellerca/provean/bin/provean/bin/provean.sh

# To edit any of the configuration
vi /data/${USER}/provean/bin/provean.sh

# /fdb/blastdb/nr
# /fdb/blastdb/v4/nr
# /fdb/blastdb/Nr_cluster/nr_cluster_seq

nohup bash /data/${USER}/provean/bin/provean.sh \
    --query ../AcrIIa4.fa \
    --variation hgvs_mutants.txt \
    --save_supporting_set AcrIIa4.sss &

module load blast
module load cd-hit

nohup bash /data/${USER}/provean/bin/provean.sh \
    --query ../AcrIIa4.fa \
    --variation hgvs_mutants.txt \
    --save_supporting_set AcrIIa4-vr.sss &


bash /data/${USER}/provean/bin/provean.sh \
    --query ../AcrIIa4.fa \
    --variation hgvs-muts.txt \
    --supporting_set AcrIIa4-vr.sss > provean-out.txt

wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ncbi-blast-2.4.0+-x64-linux.tar.gz

tar -zxvf ncbi-blast-2.4.0+-x64-linux.tar.gz

# Try with ftp version of v4
mkdir db 
cd db
wget --recursive -e robots=off --reject "indexl.html" --no-host-directories --cut-dirs=6 https://ftp.ncbi.nlm.nih.gov/blast/db/v4/ -A "nr_v4*"

```