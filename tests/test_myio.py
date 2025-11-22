import pytest
from pathlib import Path
from scripts.myio import parse_structure_pdb, parse_structure_mmcif

def test_parse_structure_pdb_extension():
    """Test that parse_structure_pdb raises ValueError for invalid extensions."""
    with pytest.raises(ValueError, match="Expected .pdb file"):
        parse_structure_pdb(Path("test.txt"))

def test_parse_structure_mmcif_extension():
    """Test that parse_structure_mmcif raises ValueError for invalid extensions."""
    with pytest.raises(ValueError, match="Expected .cif or .mmcif file"):
        parse_structure_mmcif(Path("test.txt"))
