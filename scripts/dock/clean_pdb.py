#!/usr/bin/env python3
import sys

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa


class ProteinOnly(Select):
    def accept_residue(self, residue):
        # ASSUMPTION: Removing all non-amino acid residues (including structural waters).
        # WHY IT'S UNREALISTIC GENERALLY: For most HIV-1 PR inhibitors (like Saquinavir), a conserved
        # flap water is critical for binding.
        # WHY IT'S ACCEPTABLE HERE: Darunavir is specifically designed to displace this water molecule,
        # so stripping it before docking is scientifically valid in this exact context.
        # Reference: https://pubs.acs.org/doi/10.1021/jm060561m
        return is_aa(residue, standard=False)

    def accept_atom(self, atom):
        return atom.get_altloc() == " "


input_pdb, output_pdb = sys.argv[1], sys.argv[2]

structure = PDBParser(QUIET=True).get_structure("receptor", input_pdb)

for atom in structure.get_atoms():
    if atom.is_disordered():
        # ASSUMPTION: The alternate location with the highest occupancy is the most representative conformation.
        # WHY IT IS UNREALISTIC SOMETIMES: Ligand binding can induce or prefer a lower-occupancy state. 
        # WHY IT IS ACCEPTABLE HERE: For rigid docking, we must select a single state, and highest occupancy is the standard unbiased choice.
        # Reference: https://pubs.acs.org/doi/10.1021/acs.jcim.3c00100 ; https://documentation.samson-connect.net/tutorials/prepare-protein/prepare-protein/
        highest_occ_child = max(atom.child_dict.values(), key=lambda x: x.get_occupancy())
        atom.disordered_select(highest_occ_child.get_altloc())
        atom.selected_child.set_altloc(" ")

io = PDBIO()
io.set_structure(structure)
io.save(output_pdb, ProteinOnly())
