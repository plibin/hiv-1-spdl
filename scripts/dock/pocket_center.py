from Bio.PDB import PDBParser
import numpy as np
import sys

pdb = sys.argv[1]
anchors = sys.argv[2:]  # e.g. A:25:OD1,OD2 B:25:OD1,OD2

s = PDBParser(QUIET=True).get_structure("x", pdb)
coords = []

for a in anchors:
    chain_id, resnum, atom_names = a.split(":")
    atom_names = atom_names.split(",")

    r = s[0][chain_id][(" ", int(resnum), " ")]
    for atom_name in atom_names:
        coords.append(r[atom_name].coord)

print(*np.mean(coords, axis=0))
