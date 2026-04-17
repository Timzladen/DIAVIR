#!/bin/bash

#run the assembly on CLC Genomics QCed and trimmed fastq files
conda activate spades
bash 01_assembly.sh
conda deactivate

#predict ORFs
conda activate transdecoder
TransDecoder.LongOrfs -t merged_transcripts.fasta
TransDecoder.Predict -t merged_transcripts.fasta
conda deactivate
#cluster contigs
conda activate cd-hit
bash 02_cluster.sh
conda deactivate
#run salmon read quantification
conda activate salmon
bash 03_salmon.sh
conda deactivate

#annotate contigs with eggnog
emapper.py \
  -i transcripts.fasta \
  -o eggnog_annotations \
  --cpu 64 \
  -m diamond \
  --itype CDS \
  --data_dir eggnogDB
