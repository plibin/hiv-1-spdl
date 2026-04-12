#!/bin/bash

data_dir=../data/

python query_fasta.py --base-path $data_dir --protein RT > $data_dir/RT/refs/query.fasta 
#there are certain PDBs for which there is no RESSEQ in the PDB, we provide an override for these
python query_fasta.py --base-path $data_dir --protein PR --overrides_fasta $data_dir/PR/refs/overrides.fasta > $data_dir/PR/refs/query.fasta 
python query_fasta.py --base-path $data_dir --protein IN --overrides_fasta $data_dir/IN/refs/overrides.fasta > $data_dir/IN/refs/query.fasta 

python pdb_fasta.py --base-path $data_dir --protein RT > $data_dir/RT/refs/pdb.fasta 
python pdb_fasta.py --base-path $data_dir --protein PR > $data_dir/PR/refs/pdb.fasta 
python pdb_fasta.py --base-path $data_dir --protein IN > $data_dir/IN/refs/pdb.fasta

#align
cat "$data_dir/RT/refs/hxb2.fasta" "$data_dir/RT/refs/query.fasta" "$data_dir/RT/refs/pdb.fasta" | mafft --auto - | seqkit sort > "$data_dir/RT/refs/alignment.fasta" 
cat "$data_dir/PR/refs/hxb2.fasta" "$data_dir/PR/refs/query.fasta" "$data_dir/PR/refs/pdb.fasta" | mafft --auto - | seqkit sort > "$data_dir/PR/refs/alignment.fasta" 
cat "$data_dir/IN/refs/hxb2.fasta" "$data_dir/IN/refs/query.fasta" "$data_dir/IN/refs/pdb.fasta" | mafft --auto - | seqkit sort > "$data_dir/IN/refs/alignment.fasta" 

python format_fasta.py --fasta "$data_dir/RT/refs/alignment.fasta"
python format_fasta.py --fasta "$data_dir/PR/refs/alignment.fasta"
python format_fasta.py --fasta "$data_dir/IN/refs/alignment.fasta"

#IN alignment has quite a bit of loose ends at the left side: trim the left side by replacing all these loose ends with gaps
python trim_alignment_left.py "$data_dir/IN/refs/alignment.fasta" 55 > "$data_dir/IN/refs/alignment.fasta.clean"
mv "$data_dir/IN/refs/alignment.fasta.clean" "$data_dir/IN/refs/alignment.fasta" 
