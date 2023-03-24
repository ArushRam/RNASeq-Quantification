kallisto index -i kallisto/index data/chr11_transcriptome.fasta
kallisto quant -i kallisto/index -o kallisto/output -l 200 -s 10 --single data/reads.fastq