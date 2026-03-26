import copy

import numpy as np

from scripts.plddt import per_residue_plddt
from scripts.rmsd import global_rmsd, per_residue_rmsd


def test_per_residue_rmsd_self_alignment_is_zero(sample_chain):
    pred_chain = copy.deepcopy(sample_chain)

    result = per_residue_rmsd("toy", "AGS", "AGS", sample_chain, pred_chain)

    assert result.keys() == {0, 1, 2}
    assert all(np.isclose(value, 0.0) for value in result.values())


def test_global_rmsd_self_alignment_is_zero(sample_chain):
    pred_chain = copy.deepcopy(sample_chain)

    assert np.isclose(global_rmsd("toy", "AGS", "AGS", sample_chain, pred_chain), 0.0)


def test_per_residue_plddt_reads_ca_bfactors(make_chain):
    ref_chain = make_chain("AGS")
    pred_chain = make_chain("AGS", ca_bfactors=[91.5, 82.0, 70.25])

    result = per_residue_plddt("toy", "AGS", "AGS", ref_chain, pred_chain)

    assert result == {0: 91.5, 1: 82.0, 2: 70.25}
