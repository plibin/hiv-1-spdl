import pytest
from pathlib import Path
from scripts.myio import load_preds, load_refs, parse_structure_pdb, parse_structure_mmcif


def test_parse_structure_pdb_extension():
    """Test that parse_structure_pdb raises ValueError for invalid extensions."""
    with pytest.raises(ValueError, match="Expected .pdb file"):
        parse_structure_pdb(Path("test.txt"))


def test_parse_structure_mmcif_extension():
    """Test that parse_structure_mmcif raises ValueError for invalid extensions."""
    with pytest.raises(ValueError, match="Expected .cif or .mmcif file"):
        parse_structure_mmcif(Path("test.txt"))


def test_load_refs_uses_uppercase_ids(data_root):
    refs = load_refs(data_root, "PR")

    assert "7LEG" in refs
    assert "7leg" not in refs


def test_load_preds_returns_known_prediction(data_root):
    refs = {"7LEG": parse_structure_pdb(data_root / "PR" / "refs" / "7leg.pdb")}

    preds = load_preds(refs, data_root, "PR", "AlphaFold2")

    assert set(preds) == {"7LEG"}
