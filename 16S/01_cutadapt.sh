#!/bin/bash
cd D5C2
cutadapt -g AACGTGAT -G AACCGAGA -a GTGYCAGCMGCCGCGGTAA -A CCGYCAATTYMTTTRAGTTT -o D5C2.raw_1_cutad.fastq.gz -p D5C2.raw_2_cutad.fastq.gz D5C2.raw_1.fastq.gz D5C2.raw_2.fastq.gz
cd ../D5C3
cutadapt -g AACGTGAT -G AACGCTTA -a GTGYCAGCMGCCGCGGTAA -A CCGYCAATTYMTTTRAGTTT -o D5C3.raw_1_cutad.fastq.gz -p D5C3.raw_2_cutad.fastq.gz D5C3.raw_1.fastq.gz D5C3.raw_2.fastq.gz
cd ../D5V1
cutadapt -g AAACATCG -G AGTACAAG -a GTGYCAGCMGCCGCGGTAA -A CCGYCAATTYMTTTRAGTTT -o D5V1.raw_1_cutad.fastq.gz -p D5V1.raw_2_cutad.fastq.gz D5V1.raw_1.fastq.gz D5V1.raw_2.fastq.gz
cd ../D5V2
cutadapt -g AAACATCG -G AACAACCA -a GTGYCAGCMGCCGCGGTAA -A CCGYCAATTYMTTTRAGTTT -o D5V2.raw_1_cutad.fastq.gz -p D5V2.raw_2_cutad.fastq.gz D5V2.raw_1.fastq.gz D5V2.raw_2.fastq.gz
cd ../D5V3
cutadapt -g AAACATCG -G AACCGAGA -a GTGYCAGCMGCCGCGGTAA -A CCGYCAATTYMTTTRAGTTT -o D5V3.raw_1_cutad.fastq.gz -p D5V3.raw_2_cutad.fastq.gz D5V3.raw_1.fastq.gz D5V3.raw_2.fastq.gz
cd ../D9C1
cutadapt -g AAACATCG -G AACGCTTA -a GTGYCAGCMGCCGCGGTAA -A CCGYCAATTYMTTTRAGTTT -o D9C1.raw_1_cutad.fastq.gz -p D9C1.raw_2_cutad.fastq.gz D9C1.raw_1.fastq.gz D9C1.raw_2.fastq.gz
cd ../D9C2
cutadapt -g ATGCCTAA -G AGTACAAG -a GTGYCAGCMGCCGCGGTAA -A CCGYCAATTYMTTTRAGTTT -o D9C2.raw_1_cutad.fastq.gz -p D9C2.raw_2_cutad.fastq.gz D9C2.raw_1.fastq.gz D9C2.raw_2.fastq.gz
cd ../D9C3
cutadapt -g ATGCCTAA -G AACAACCA -a GTGYCAGCMGCCGCGGTAA -A CCGYCAATTYMTTTRAGTTT -o D9C3.raw_1_cutad.fastq.gz -p D9C3.raw_2_cutad.fastq.gz D9C3.raw_1.fastq.gz D9C3.raw_2.fastq.gz
cd ../D9V1
cutadapt -g ATGCCTAA -G AACCGAGA -a GTGYCAGCMGCCGCGGTAA -A CCGYCAATTYMTTTRAGTTT -o D9V1.raw_1_cutad.fastq.gz -p D9V1.raw_2_cutad.fastq.gz D9V1.raw_1.fastq.gz D9V1.raw_2.fastq.gz
cd ../D9V2
cutadapt -g ATGCCTAA -G AACGCTTA -a GTGYCAGCMGCCGCGGTAA -A CCGYCAATTYMTTTRAGTTT -o D9V2.raw_1_cutad.fastq.gz -p D9V2.raw_2_cutad.fastq.gz D9V2.raw_1.fastq.gz D9V2.raw_2.fastq.gz
cd ../D9V3
cutadapt -g AGTGGTCA -G AGTACAAG -a GTGYCAGCMGCCGCGGTAA -A CCGYCAATTYMTTTRAGTTT -o D9V3.raw_1_cutad.fastq.gz -p D9V3.raw_2_cutad.fastq.gz D9V3.raw_1.fastq.gz D9V3.raw_2.fastq.gz
cd ../Int
cutadapt -g AACGTGAT -G AGTACAAG -a GTGYCAGCMGCCGCGGTAA -A CCGYCAATTYMTTTRAGTTT -o Int.raw_1_cutad.fastq.gz -p Int.raw_2_cutad.fastq.gz Int.raw_1.fastq.gz Int.raw_2.fastq.gz
