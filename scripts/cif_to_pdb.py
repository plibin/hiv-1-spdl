import sys
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBIO

def cif_to_pdb(cif_path: Path):
    # Parse the CIF structure
    parser = MMCIFParser(QUIET=True)
    structure_id = cif_path.stem
    structure = parser.get_structure(structure_id, cif_path)

    # Output filename (same name, .pdb extension)
    pdb_path = cif_path.with_suffix(".pdb")

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(pdb_path))

def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: python cif_to_pdb.py <path/to/file.cif>\n")
        sys.exit(1)

    cif_to_pdb(Path(sys.argv[1]))

if __name__ == "__main__":
    main()
