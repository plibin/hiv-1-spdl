#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument("--fasta", required=True)
parser.add_argument("--seqname", required=True)
parser.add_argument("--pos", type=int, required=True)  # 1-based
args = parser.parse_args()

records = list(SeqIO.parse(args.fasta, "fasta"))

for rec in records:
    if rec.id == args.seqname:
        seq = str(rec.seq)
        print(seq)
        rec.seq = Seq(seq[:args.pos - 1] + "-" + seq[args.pos:])
        print(rec.seq)
        break
else:
    raise ValueError(f"Sequence {args.seqname} not found")

SeqIO.write(records, args.fasta, "fasta-2line")
