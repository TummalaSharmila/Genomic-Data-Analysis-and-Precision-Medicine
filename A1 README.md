## **COURSE NAME AND CODE:** Genomic Data Analysis and Precision Medicine
### **Programmer:** Sharmila Tummala

Date: 10/17/2025

### **Programming Language:** 
Unix, Python


### **Description:** 
This approach uses short-read Illumina sequencing  genome data of Escherichia coli.  After the data was preprocessed using  FastQC  and Trim Galore  for quality control and trimming, it was assembled using two de novo genome assembly methods, Velvet and Oases, with different k-mer sizes (31, 41, 51, 61, and 71).  The goal was to evaluate the performance of both technologies and determine the optimal genome assembly parameters.

### TOOLS USED
SRA Toolkit: To download sequencing data from the NCBI SRA database.
Trim Galore: For trimming adapter sequences and low-quality bases.
FastQC: For evaluating the quality of the raw and trimmed reads.
Velvet: A de novo genomic assembler optimized for short-read data.
Oases: A transcriptome assembler based on Velvet but applied here for genome assembly.
QUAST: A tool used for evaluating the quality of genome assemblies.



### Packages needed

