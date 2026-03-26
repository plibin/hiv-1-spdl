#!/bin/bash

data_dir=../data/

python query_fasta.py --base-path $data_dir --protein RT > $data_dir/RT/refs/query.fasta 
python query_fasta.py --base-path $data_dir --protein PR > $data_dir/PR/refs/query.fasta 
python query_fasta.py --base-path $data_dir --protein IN > $data_dir/IN/refs/query.fasta 

python pdb_fasta.py --base-path $data_dir --protein RT > $data_dir/RT/refs/pdb.fasta 
python pdb_fasta.py --base-path $data_dir --protein PR > $data_dir/PR/refs/pdb.fasta 
python pdb_fasta.py --base-path $data_dir --protein IN > $data_dir/IN/refs/pdb.fasta

#align
cat "$data_dir/RT/refs/query.fasta" "$data_dir/RT/refs/pdb.fasta" | mafft --auto - | seqkit sort > "$data_dir/RT/refs/alignment.fasta" 
cat "$data_dir/PR/refs/query.fasta" "$data_dir/PR/refs/pdb.fasta" | mafft --auto - | seqkit sort > "$data_dir/PR/refs/alignment.fasta" 
cat "$data_dir/IN/refs/query.fasta" "$data_dir/IN/refs/pdb.fasta" | mafft --auto - | seqkit sort > "$data_dir/IN/refs/alignment.fasta" 
