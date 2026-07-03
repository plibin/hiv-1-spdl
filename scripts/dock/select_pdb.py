from sys import argv
from pathlib import Path
from shutil import copy2
from Bio.PDB import PDBParser, PDBIO, Select

pdb, outdir = Path(argv[1]), Path(argv[2])
outdir.mkdir(exist_ok=True)


class KeepAB(Select):
    def accept_chain(self, chain):
        return chain.id in {"A", "B"}

def first_standard_residue_number(chain):
    for residue in chain:
        if residue.id[0] == " ":
            return residue.id[1]
    return None

def renumber_chain_b(model):
    i = 1
    for residue in model["B"]:
        if residue.id[0] == " ":
            residue.id = (" ", i, " ")
            i += 1

s = PDBParser(QUIET=True).get_structure(pdb.stem, pdb)
chains = {c.id for c in s[0]}

if 1 <= len(chains) <= 99 and {"A", "B"} <= chains:
    out = outdir / pdb.name
    
    #for ESMFold, the PDBs have numbering that does not start at 1, 
    #for the B-chain -> fix it
    if first_standard_residue_number(s[0]["B"]) != 1:
      renumber_chain_b(s[0])

    io = PDBIO()
    io.set_structure(s)
    io.save(str(out), KeepAB())
