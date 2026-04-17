#!/usr/bin/env bash
set -euo pipefail

# ==========================================
# Build MMETSP Pseudo-nitzschia reference,
# annotate, rename eggNOG IDs, and quantify
# with Salmon
# ==========================================

# 1. Prepare combined MMETSP reference
cd /DATA/scratch/timtd/transcriptomes/MMETSP/fasta

# Combine FASTAs
cat *.fasta > combined_pn_transcriptome.fasta

# Dereplicate similar transcripts
conda activate cdhit
cd-hit-est \
  -i combined_pn_transcriptome.fasta \
  -o cleaned_pn_transcriptome.fasta \
  -c 0.95 \
  -n 10

# 2. Annotate transcriptome
conda activate dammit
dammit annotate cleaned_pn_transcriptome.fasta \
  --busco-group stramenopiles \
  --n_threads 32



#!/bin/bash

# Set paths
THREADS=64
DAMMIT_DIR="cleaned_pn_transcriptome.fasta.dammit"
ORF_FASTA="$DAMMIT_DIR/cleaned_pn_transcriptome.fasta.transdecoder_dir/longest_orfs.pep"
PFAM_DB="$HOME/dammit_db/Pfam-A.hmm"
OUTPUT="$DAMMIT_DIR/cleaned_pn_transcriptome.fasta.transdecoder_dir/longest_orfs.pep.pfam.domtblout"

# Run hmmscan in parallel
echo "Running hmmscan with $THREADS threads..."
hmmscan --cpu $THREADS \
        --domtblout "$OUTPUT" \
        -E 1e-5 \
        "$PFAM_DB" \
        "$ORF_FASTA"

# Confirm result
echo "hmmscan complete. Results saved to:"
echo "$OUTPUT"



# 3. Run eggNOG on the same cleaned reference
conda activate eggnog
emapper.py \
  -i cleaned_pn_transcriptome.fasta \
  -o pn_eggnog \
  --cpu 64 \
  -m diamond

#dammit renames the transcripts so the eggnog annotation needs to be renamed as well
python3 rename-eggnog.py

# 5. Build Salmon index on the cleaned transcriptome
conda activate salmon
salmon index \
  -t cleaned_pn_transcriptome.fasta \
  -i salmon_index_cleaned_pn \
  -k 31

# 6. Quantify samples
READS_DIR="/DATA/scratch/timtd/trimmed"
OUT_DIR="/DATA/scratch/timtd/transcriptomes/MMETSP/fasta/salmon_quant_renamed"
mkdir -p "${OUT_DIR}"

for sample in DIAVIR1 DIAVIR2 DIAVIR3 DIAVIR4 DIAVIR5 DIAVIR6 DIAVIR7 DIAVIR8 DIAVIR9 DIAVIR10 DIAVIR11 DIAVIR12 DIAVIR13
do
  salmon quant \
    -i salmon_index_cleaned_pn \
    -l A \
    -1 "${READS_DIR}/${sample}_R1.fastq.gz" \
    -2 "${READS_DIR}/${sample}_R2.fastq.gz" \
    -p 16 \
    --validateMappings \
    -o "${OUT_DIR}/${sample}"
done

