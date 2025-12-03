# create a txt file with all SRA accession numbers
cat > SraAccList.txt <<EOF
SRR23108034
SRR23108035
SRR23108036
SRR23108037
SRR23108038
SRR23108039
SRR23108040
SRR23108041
SRR23108042
SRR23108043
EOF

## Prefetch
module load sra-toolkit
prefetch --option-file SraAccList.txt -O SRR/

## convert to fastq
for acc in $(cat SraAccList.txt)
do
  fasterq-dump $acc --threads 8 --outdir /N/scratch/shatumma/metatx_project/fastq
done

## Compress FASTQ files
cd /N/scratch/shatumma/metatx_project/fastq
gzip *.fastq

## run QC check
module load fastqc
mkdir -p qc

fastqc fastq/*.gz -o qc

