#!/usr/bin/env python3

import sys
from pathlib import Path

if len(sys.argv) != 4:
    sys.exit("Usage: rename_fasta_and_make_namemap.py <input.fasta> <output.renamed.fasta> <output.namemap.csv>")

in_fasta = Path(sys.argv[1])
out_fasta = Path(sys.argv[2])
out_csv = Path(sys.argv[3])

counter = 0

with in_fasta.open() as fin, out_fasta.open("w") as fout, out_csv.open("w") as fmap:
    fmap.write("original_id,renamed_id\n")

    for line in fin:
        if line.startswith(">"):
            original_header = line[1:].strip()
            original_id = original_header.split()[0]
            renamed_id = f"Transcript_{counter}"
            counter += 1

            fout.write(f">{renamed_id}\n")
            fmap.write(f"{original_id},{renamed_id}\n")
        else:
            fout.write(line)

print(f"Renamed {counter} transcripts.")
print(f"Renamed FASTA: {out_fasta}")
print(f"Name map: {out_csv}")
