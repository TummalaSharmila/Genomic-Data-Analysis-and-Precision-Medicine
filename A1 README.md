## **COURSE NAME AND CODE:** Genomic Data Analysis and Precision Medicine
### **Programmer:** Sharmila Tummala

Date: 10/17/2025

### **Programming Language:** 
Unix, Python


### **Description:** 
This approach uses short-read Illumina sequencing  genome data of Escherichia coli.  After the data was preprocessed using  FastQC  and Trim Galore  for quality control and trimming, it was assembled using two de novo genome assembly methods, Velvet and Oases, with different k-mer sizes (31, 41, 51, 61, and 71).  The goal was to evaluate the performance of both technologies and determine the optimal genome assembly parameters.

### TOOLS USED
SRA Toolkit: To download sequencing data from the NCBI SRA database.
FastQC: For evaluating the quality of the raw and trimmed reads.
Trim Galore: For trimming adapter sequences and low-quality bases.
Velvet: A de novo genomic assembler optimized for short-read data.
Oases: A transcriptome assembler based on Velvet but applied here for genome assembly.
QUAST: A tool used for evaluating the quality of genome assemblies.
## Workflow
I created an envirnoment named Assignment1 to do this task using conda in HPC by loading conda as follows
module load conda #loads conda
conda create --name Assignment1 #creates environment with the name Assignment1
conda install bioconda::velvet oases
pip install quast # installs QUAST v5.2.0
mkdir /..path/Assignment1
cd /..path/Assignment1

### Data Extraction
module load sra-toolkit #loads package
fasterq-dump --split-files SRR21904868 #downloads the 2 sequencing files

### QC control
module load fastqc
mkdir fastqc
fastqc -o fastqc SRR30659984_1.fastq SRR30659984_2.fastq

### trimming
module load python
pip install cutadapt
module load trimgalore
module load trimgalore
trim_galore --phred33 --fastqc --paired SRR21904868_1.fastq SRR21904868_2.fastq -o trimming

mkdir -p assemblies

### Velvet Assembly
to run velvet on the chosen k-mers run the below code in command line (The lopp run the velvet for each k-mer and stores them in respective folders) 

for k in 31 41 51 61 71;
do   echo "Running Velvet with k=$k";      mkdir -p assemblies/velvet_k$k;      
velveth assemblies/velvet_k$k $k -shortPaired -fastq -separate SRR21904868_1.fastq SRR21904868_2.fastq;      
velvetg assemblies/velvet_k$k -exp_cov auto -cov_cutoff auto;      
echo "Velvet assembly complete for k=$k"; 
done

This generated the following files for each k-mer:
contigs.fa
Graph2
LastGraph
Log
PreGraph
Roadmaps
Sequences
stats.txt

### Oases Assembly
to run oases on the chosen k-mers run the below code in command line (The lopp run the oases for each k-mer and stores them in respective folders) 

for k in 31 41 51 61 71; do   
echo "Running Oases with k-mer = $k";   
mkdir -p assemblies/oases_k$k;   
cp -r assemblies/velvet_k$k/* assemblies/oases_k$k/;   
cd assemblies/oases_k$k;   
oases .;   
cd ../../; 
done

This generated the following files for each k-mer:
contigs-orderng.txt
contigs.fa
Graph2
LastGraph
Log
PreGraph
Roadmaps
Sequences
stats.txt
transcripts.fa

Now we have all the k-mer assemblies from velvet and oases. I am using quast to analyse which assembly tool is best by comparing various aspects

### Quast Analysis
mkdir -p quast_results

#quast for velvet assemblies

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

### Conclusion:
Best Tool: Oases outperformed Velvet by a large margin in terms of contig length, N50, and total assembled length, reflecting its suitability for transcriptome like assemblies.

Optimal k-mer Size: k-mer 51 was optimal for Oases, offering the best balance of contiguity and completeness, while Velvet failed to achieve quality assemblies across tested k-mer sizes.

