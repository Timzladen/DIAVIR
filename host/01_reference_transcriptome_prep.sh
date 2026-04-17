#!/bin/bash

cd /DATA/scratch/timtd/transcriptomes/MMETSP/fasta

# Combine FASTAs
cat *.fasta > combined_pn_transcriptome.fasta

conda activate cdhit
#dereplicate similar transcripts
cd-hit-est -i combined_pn_transcriptome.fasta -o cleaned_pn_transcriptome.fasta -c 0.95 -n 10

