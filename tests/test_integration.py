import pytest
from Bio import SeqIO

from scripts.myio import _first_chain, _first_model, parse_structure_pdb
from scripts.plddt import per_residue_plddt
from scripts.rmsd import global_rmsd
from scripts.tm import global_tm
from scripts.rmsd import per_residue_rmsd


@pytest.fixture(scope="session")
def pr_7leg_alignment(data_root):
    records = SeqIO.to_dict(SeqIO.parse(data_root / "PR" / "refs" / "alignment.fasta", "fasta"))
    return str(records["7LEG_pdb"].seq), str(records["7LEG"].seq)


@pytest.fixture(scope="session")
def pr_7leg_chains(data_root):
    ref_structure = parse_structure_pdb(data_root / "PR" / "refs" / "7leg.pdb")
    pred_structure = parse_structure_pdb(data_root / "PR" / "AlphaFold2" / "7leg.pdb")
    return _first_chain(_first_model(ref_structure)), _first_chain(_first_model(pred_structure))


def test_alignment_pair_lengths_match(pr_7leg_alignment):
    ref_align, pred_align = pr_7leg_alignment
    assert len(ref_align) == len(pred_align)
    assert len(ref_align) > 0


def test_rmsd_integration(pr_7leg_alignment, pr_7leg_chains):
    ref_align, pred_align = pr_7leg_alignment
    ref_chain, pred_chain = pr_7leg_chains

    rmsd_by_pos = per_residue_rmsd("7LEG", ref_align, pred_align, ref_chain, pred_chain)

    assert rmsd_by_pos
    assert all(pos >= 0 for pos in rmsd_by_pos)
    assert all(value >= 0.0 for value in rmsd_by_pos.values())


def test_global_rmsd_integration(pr_7leg_alignment, pr_7leg_chains):
    ref_align, pred_align = pr_7leg_alignment
    ref_chain, pred_chain = pr_7leg_chains

    rmsd_value = global_rmsd("7LEG", ref_align, pred_align, ref_chain, pred_chain)

    assert rmsd_value >= 0.0


def test_plddt_integration(pr_7leg_alignment, pr_7leg_chains):
    ref_align, pred_align = pr_7leg_alignment
    ref_chain, pred_chain = pr_7leg_chains

    plddt_scores = per_residue_plddt("7LEG", ref_align, pred_align, ref_chain, pred_chain)

    assert plddt_scores
    assert all(0.0 <= score <= 100.0 for score in plddt_scores.values())


def test_tm_integration(pr_7leg_alignment, pr_7leg_chains):
    ref_align, pred_align = pr_7leg_alignment
    ref_chain, pred_chain = pr_7leg_chains

    tm_score = global_tm("7LEG", ref_align, pred_align, ref_chain, pred_chain)

    assert 0.0 <= tm_score <= 1.0
