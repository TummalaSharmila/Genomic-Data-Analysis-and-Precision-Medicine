###create environment
conda create --name Assignment1
conda activate Assignment1
conda install -c bioconda:: velvet oasis
mkdir Assignment1
cd Assignment1

###download the SRR21904868 files
fasterq-dump --split-files SRR21904868 

###Quality control 
mkdir fastqc
fastqc -o fastqc SRR30659984_1.fastq SRR30659984_2.fastq

#trimming 
pip install cutadapt
module load trimgalore
trim_galore --phred33 --fastqc --paired SRR21904868_1.fastq SRR21904868_2.fastq -o trimming

mkdir -p assemblies

#vevet
for k in 31 41 51 61 71;
do   echo "Running Velvet with k=$k";      mkdir -p assemblies/velvet_k$k;      
velveth assemblies/velvet_k$k $k -shortPaired -fastq -separate SRR21904868_1.fastq SRR21904868_2.fastq;      
velvetg assemblies/velvet_k$k -exp_cov auto -cov_cutoff auto;      
echo "Velvet assembly complete for k=$k"; 
done

#oasis
for k in 31 41 51 61 71; do   
echo "Running Oases with k-mer = $k";   
mkdir -p assemblies/oases_k$k;   
cp -r assemblies/velvet_k$k/* assemblies/oases_k$k/;   
cd assemblies/oases_k$k;   
oases .;   
cd ../../; 
done

#quast for velvet assemblies
mkdir -p quast_results

for k in 31 41 51 61 71
do
    echo "Running QUAST for Velvet k=$k..."
    quast.py assemblies/velvet_k${k}/contigs.fa \
        -o quast_results/velvet_k${k} \
        --min-contig 200 --threads 8
done

#quast for oases assemblies

for k in 31 41 51 61 71
do
    quast.py assemblies/oases_k${k}/transcripts.fa -o quast_results/oases_k${k} --min-contig 200 --threads 8
done
