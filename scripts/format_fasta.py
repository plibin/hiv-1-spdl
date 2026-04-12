#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument("--fasta", required=True)
args = parser.parse_args()

records = list(SeqIO.parse(args.fasta, "fasta"))

SeqIO.write(records, args.fasta, "fasta-2line")
