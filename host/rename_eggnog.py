#!/usr/bin/env python3

import os
import csv
import pandas as pd

# Input files
NAMEMAP_PATH = '/users/timtd/transcriptomes/MMETSP/fasta/cleaned_pn_transcriptome.fasta.dammit/cleaned_pn_transcriptome.fasta.dammit.namemap.csv'
EGGNOG_PATH = 'pn_eggnog.emapper.annotations'
OUTPUT_PATH = 'pn_eggnog.renamed.tsv'

# Load ID mapping
id_map = {}
with open(NAMEMAP_PATH, newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if row and not row[0].startswith('#') and not row[0].startswith('original'):
            original_id, renamed_id = row
            clean_id = original_id.strip().split()[0]
            id_map[clean_id] = renamed_id.strip()

# Detect header line in eggNOG file
with open(EGGNOG_PATH, 'r') as f:
    lines = f.readlines()

# Get column names from 5th line (index 4)
header = lines[4].strip().split('\t')

# Read the data starting from line 5
data = pd.read_csv(EGGNOG_PATH, sep='\t', comment='#', header=None, skiprows=5, names=header)

# Rename the query column
data['#query'] = data['#query'].apply(lambda qid: id_map.get(qid.strip(), qid.strip()))

# Save output
data.to_csv(OUTPUT_PATH, sep='\t', index=False)
print(f" ^|^e Renamed eggNOG annotations saved to: {OUTPUT_PATH}")
