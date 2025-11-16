# Genomic-Analysis-Assignment-2
## Author: Sharmila Tummala
# Title: Differential Expression Analysis of SARS-CoV-2 and Mock Infected Cells Across Time Points
The goal of this assignment is to compare the gene expression of human respiratory cells infected with SARS-CoV-2 with mock (control) samples at two different time points (24 hours and 72 hours). Data collection and preprocessing, alignment, differential expression and enrichment analysis, and sample and condition annotation of expression data are the primary phases.

# Downloading data

module load sra-toolkit   #v3.0.5

mkdir -p SRR data/fastq

## loading all the given data file's accession numbers in a .txt file
cat > SraAccList.txt <<EOF
SRR22269883
SRR22269882
SRR22269881
SRR22269880
SRR22269879
SRR22269878
SRR22269877
SRR22269876
SRR22269875
SRR22269874
SRR22269873
SRR22269872
EOF

# Download SRA files
prefetch --option-file SraAccList.txt -O SRR/

# Convert to fastq
for acc in $(cat SraAccList.txt); do
  fasterq-dump $acc -O data/fastq --threads 8
done


# Quality Control
module load fastqc   #v0.12.1

mkdir -p qc
fastqc data/fastq/*.fastq -o qc

pip install cutadapt #v4.9
module load trimgalore #v0.6.10

mkdir -p trimming
for fq in data/fastq/*.fastq; do
  trim_galore -q 30  --length 16 --phred33 --fastqc "$fq" -o trimming
done


# Reference Setup

mkdir -p reference
cd reference
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz   #release 115 latest one from ensembl
gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz

cd ..

# Indexing

module load star # v2.7.11a
mkdir -p hg38_115
#buidling index
STAR --runMode genomeGenerate \
     --genomeDir hg38_115 \
     --genomeFastaFiles reference/Homo_sapiens.GRCh38.dna.toplevel.fa \
     --runThreadN 4

# Alignment

mkdir -p alignment
SRA_accessions=("SRR22269883" "SRR22269882" "SRR22269881" "SRR22269880" "SRR22269879" "SRR22269878" "SRR22269877" "SRR22269876" "SRR22269875" "SRR22269874" "SRR22269873" "SRR22269872")

STAR --genomeDir hg38_115 \
     --readFilesIn ${SRA_Accession}.fastq \
     --outFileNamePrefix alignment/"${SRA_Accession}." \
     --outSAMtype BAM SortedByCoordinate \
     --runThreadN 4

# Quantification
## Downloading .gtf file
cd reference 
wget https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.115.chr.gtf.gz #release 115
gunzip Homo_sapiens.GRCh38.115.chr.gtf.gz

cd ..

module load subread
### run the below command for all the bam files to get counts.txt
featureCounts -a Homo_sapiens.GRCh38.115.chr.gtf -o gene_counts.txt alignments/*.sorted.bam

# Extract gene counts of each SRA_Accession number and make a gene_counts.csv file 
Manual combination of all txt in one csv in excel
>open excel >data >get data >select the .txt files >copy the counts for each SRA_Accession into a new csv with geneID column and get a complete gene_counts.csv containing all the SRA_Accession counts for the gene IDs

# Running R Script to Differential Expression Analysis of SARS-CoV-2 and Mock Infected Cells Across Time Points
I am uploading the .R file indicating all the steps done there


# Interpretation of Differential Expression Results
