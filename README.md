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
->open excel ->data ->get data ->select the .txt files ->copy the counts for each SRA_Accession into a new csv with geneID column and get a complete gene_counts.csv containing all the SRA_Accession counts for the gene IDs
submitting the **gene_counts.csv** with all the gene IDs and there counts in across each condition/sample in each of these two sets.
The file contains the expression levels for all genes, and miRNAs are included and detectable by their gene ID.

# Running R Script to Differential Expression Analysis of SARS-CoV-2 and Mock Infected Cells Across Time Points
I am uploading the .R file indicating all the steps done there


# Interpretation of Differential Expression Results
## Mock vs SARS-CoV-2
upregulated_gene_ids

 [1] "ENSG00000202031" "ENSG00000201321" "ENSG00000234741" "ENSG00000143226" "ENSG00000202400"
 [6] "ENSG00000265706" "ENSG00000202164" "ENSG00000222345" "ENSG00000199609" "ENSG00000202054"
[11] "ENSG00000221500" "ENSG00000305361" "ENSG00000201772" "ENSG00000252316" "ENSG00000200852"
[16] "ENSG00000285776" "ENSG00000221716" "ENSG00000233998" "ENSG00000238886" "ENSG00000294255"
[21] "ENSG00000222477" "ENSG00000287264" "ENSG00000255717" "ENSG00000255008" "ENSG00000202314"
[26] "ENSG00000199535" "ENSG00000221164" "ENSG00000302642" "ENSG00000207031" "ENSG00000251898"
[31] "ENSG00000259932" "ENSG00000261441" "ENSG00000199568" "ENSG00000258947" "ENSG00000275084"
[36] "ENSG00000212163" "ENSG00000238793" "ENSG00000199874" "ENSG00000238531" "ENSG00000209702"
[41] "ENSG00000201025" "ENSG00000210077" "ENSG00000210082"

downregulated_gene_ids
 [1] "ENSG00000188452" "ENSG00000283203" "ENSG00000115317" "ENSG00000244710" "ENSG00000106511"
 [6] "ENSG00000207175" "ENSG00000244642" "ENSG00000198879" "ENSG00000043591" "ENSG00000202347"
[11] "ENSG00000185352" "ENSG00000243700" "ENSG00000242894" "ENSG00000222076" "ENSG00000179776"
[16] "ENSG00000134438" "ENSG00000101003"

A total of 60 genes were significantly differentially expressed between mock and SARS-CoV-2 infected cells (all time points combined).

Of these, 43 genes were upregulated (higher expression in SARS-CoV-2) and 17 genes were downregulated (lower expression in SARS-CoV-2).

The mean log2 fold change was 1.12 and the median was 1.87, indicating that, on average, gene expression is increased in SARS-CoV-2 infected cells compared to mock.

The adjusted p-values (min: 2.8×10⁻⁶, max: 0.045) confirm the statistical significance of these changes.

## SARS-CoV-2 24H vs 72H
upregulated_gene_ids_24H_72H
 [1] "ENSG00000199047" "ENSG00000281706" "ENSG00000201943" "ENSG00000272460" "ENSG00000201143"
 [6] "ENSG00000138606" "ENSG00000283160" "ENSG00000131747" "ENSG00000141736" "ENSG00000274713"

 downregulated_gene_ids_24H_72H
 [1] "ENSG00000212283" "ENSG00000119772" "ENSG00000212452" "ENSG00000169429" "ENSG00000131711"
 [6] "ENSG00000197081" "ENSG00000244710" "ENSG00000223259" "ENSG00000207175" "ENSG00000244642"
[11] "ENSG00000238965" "ENSG00000221514"

Within SARS-CoV-2 infected cells, 22 genes were significantly differentially expressed between 24 hours and 72 hours post-infection.

Of these, 10 genes were upregulated (higher expression at 24H) and 12 genes were downregulated (lower expression at 24H).

The mean log2 fold change was 0.03 and the median was -1.41, indicating that, on average, gene expression is slightly decreased at 24H compared to 72H.

The adjusted p-values (min: 1.3×10⁻¹², max: 0.037) confirm the statistical significance of these changes.

## Biological Implication
The comparison between mock and SARS-CoV-2 infected cells reveals a strong transcriptional response to infection, with a majority of genes being upregulated.

The comparison between 24H and 72H in SARS-CoV-2 infected cells suggests a dynamic change in gene expression over time, with a slight trend toward decreased expression at the earlier time point.
These results suggest that SARS-CoV-2 infection triggers a robust host transcriptional response, characterized by the upregulation of immune and inflammatory genes, which may reflect the activation of antiviral defenses. The dynamic changes in gene expression over time indicate a shifting host response, possibly reflecting viral replication cycles and adaptation. Such transcriptional profiles could serve as potential biomarkers for infection severity or progression and may inform the development of targeted therapies or diagnostic strategies for COVID-19.
