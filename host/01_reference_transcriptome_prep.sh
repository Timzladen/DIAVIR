#!/usr/bin/env bash
set -euo pipefail

THREADS=64
WORKDIR="/DATA/scratch/timtd/transcriptomes/MMETSP/fasta"
READS_DIR="/DATA/scratch/timtd/trimmed"

cd "${WORKDIR}"

# 1. Combine MMETSP transcriptomes
cat *.fasta > combined_pn_transcriptome.fasta

# 2. Dereplicate similar transcripts
conda activate cdhit
cd-hit-est \
  -i combined_pn_transcriptome.fasta \
  -o cleaned_pn_transcriptome.fasta \
  -c 0.95 \
  -n 10

# 3. Rename transcript IDs and create namemap
python3 rename_fasta_and_make_namemap.py \
  cleaned_pn_transcriptome.fasta \
  cleaned_pn_transcriptome.renamed.fasta \
  cleaned_pn_transcriptome.namemap.csv

# 4. Run TransDecoder explicitly
conda activate transdecoder
TransDecoder.LongOrfs -t cleaned_pn_transcriptome.renamed.fasta
TransDecoder.Predict  -t cleaned_pn_transcriptome.renamed.fasta

# 5. Run hmmscan on predicted peptides
conda activate hmmer
PFAM_DB="$HOME/dammit_db/Pfam-A.hmm"
PEP_FASTA="cleaned_pn_transcriptome.renamed.fasta.transdecoder_dir/longest_orfs.pep"
PFAM_OUT="cleaned_pn_transcriptome.renamed.fasta.transdecoder_dir/longest_orfs.pep.pfam.domtblout"

hmmscan \
  --cpu "${THREADS}" \
  --domtblout "${PFAM_OUT}" \
  -E 1e-5 \
  "${PFAM_DB}" \
  "${PEP_FASTA}"

# 6. Run eggNOG on renamed transcriptome
conda activate eggnog
emapper.py \
  -i cleaned_pn_transcriptome.renamed.fasta \
  -o pn_eggnog \
  --cpu "${THREADS}" \
  -m diamond

# 7. Rename eggNOG query IDs to match renamed transcript IDs
python3 rename-eggnog.py

# 8. Build Salmon index on renamed transcriptome
conda activate salmon
salmon index \
  -t cleaned_pn_transcriptome.renamed.fasta \
  -i salmon_index_cleaned_pn_renamed \
  -k 31

# 9. Quantify all samples
OUT_DIR="${WORKDIR}/salmon_quant_renamed"
mkdir -p "${OUT_DIR}"

for sample in DIAVIR1 DIAVIR2 DIAVIR3 DIAVIR4 DIAVIR5 DIAVIR6 DIAVIR7 DIAVIR8 DIAVIR9 DIAVIR10 DIAVIR11 DIAVIR12 DIAVIR13
do
  salmon quant \
    -i salmon_index_cleaned_pn_renamed \
    -l A \
    -1 "${READS_DIR}/${sample}_R1.fastq.gz" \
    -2 "${READS_DIR}/${sample}_R2.fastq.gz" \
    -p 16 \
    --validateMappings \
    -o "${OUT_DIR}/${sample}"
done
