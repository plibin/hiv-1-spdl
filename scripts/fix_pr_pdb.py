import sys
from Bio.PDB import PDBParser, PDBIO, Select

class MonomerSelect(Select):
    def accept_residue(self, r):
        return r.get_id()[0] == ' ' and r.get_id()[1] < 1000

pdb_file = sys.argv[1]
parser = PDBParser(QUIET=True)
s = parser.get_structure('s', pdb_file)
io = PDBIO()
io.set_structure(s)
io.save(pdb_file, MonomerSelect())
