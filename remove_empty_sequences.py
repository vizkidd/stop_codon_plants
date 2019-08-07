#!/usr/bin/env python
from Bio import SeqIO
import sys

file=sys.argv[1]

for current_seq in SeqIO.parse(file, "fasta"):
    if str(current_seq.seq) != '':
        print ">", current_seq.id, "\n", str(current_seq.seq)
