#!/bin/bash

data_dir=../data/

python pdb_fasta.py --base-path $data_dir --protein RT > $data_dir/RT/refs/pdb.fasta 
python pdb_fasta.py --base-path $data_dir --protein PR > $data_dir/PR/refs/pdb.fasta 
python pdb_fasta.py --base-path $data_dir --protein IN > $data_dir/IN/refs/pdb.fasta 
