"""
Tests that verify the FASTA files produced by fasta.sh are present,
internally consistent, and correctly aligned for every protein in the dataset.

Checks performed per protein (PR, IN, RT):
  1. File existence  – query.fasta, pdb.fasta and alignment.fasta are present.
  2. Non-empty       – every file contains at least one sequence.
  3. No duplicates   – IDs are unique within each file.
  4. Cross-file ID coverage
       * Every query ID  has a matching <ID>_pdb  entry in pdb.fasta.
       * Every query ID  appears as <ID>      in alignment.fasta.
       * Every pdb ID    appears as <ID>_pdb  in alignment.fasta.
       * No extra IDs in alignment.fasta that are absent from query/pdb.
  5. Alignment pair integrity – for each (query, pdb) pair in the alignment:
       * Both rows have the same length.
       * Neither row is entirely gaps.
       * At least one column has a non-gap character in both rows.
"""

import pytest
from Bio import SeqIO

import config


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

GAP = "-"


def _read_fasta(path):
    """Return an ordered dict {id: str(seq)} from a FASTA file."""
    return {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}


def _query_ids(query_dict):
    return set(query_dict.keys())


def _pdb_ids(pdb_dict):
    """Strip the '_pdb' suffix to get bare IDs."""
    return {k.removesuffix("_pdb") for k in pdb_dict.keys()}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session", params=config.proteins())
def protein_fastas(request, data_root):
    """For each protein, return (protein, query_dict, pdb_dict, align_dict)."""
    protein = request.param
    refs_dir = data_root / protein / "refs"
    query = _read_fasta(refs_dir / "query.fasta")
    pdb   = _read_fasta(refs_dir / "pdb.fasta")
    align = _read_fasta(refs_dir / "alignment.fasta")
    return protein, query, pdb, align


# ---------------------------------------------------------------------------
# 1. File existence
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("protein", config.proteins())
def test_fasta_files_exist(protein, data_root):
    refs_dir = data_root / protein / "refs"
    for fname in ("query.fasta", "pdb.fasta", "alignment.fasta"):
        assert (refs_dir / fname).exists(), f"{protein}/{fname} is missing"


# ---------------------------------------------------------------------------
# 2. Non-empty
# ---------------------------------------------------------------------------

def test_query_fasta_non_empty(protein_fastas):
    protein, query, _, _ = protein_fastas
    assert len(query) > 0, f"{protein}: query.fasta is empty"


def test_pdb_fasta_non_empty(protein_fastas):
    protein, _, pdb, _ = protein_fastas
    assert len(pdb) > 0, f"{protein}: pdb.fasta is empty"


def test_alignment_fasta_non_empty(protein_fastas):
    protein, _, _, align = protein_fastas
    assert len(align) > 0, f"{protein}: alignment.fasta is empty"


# ---------------------------------------------------------------------------
# 3. No duplicate IDs
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("protein", config.proteins())
def test_no_duplicate_query_ids(protein, data_root):
    path = data_root / protein / "refs" / "query.fasta"
    ids = [r.id for r in SeqIO.parse(path, "fasta")]
    assert len(ids) == len(set(ids)), f"{protein}: duplicate IDs in query.fasta: {[i for i in ids if ids.count(i) > 1]}"


@pytest.mark.parametrize("protein", config.proteins())
def test_no_duplicate_pdb_ids(protein, data_root):
    path = data_root / protein / "refs" / "pdb.fasta"
    ids = [r.id for r in SeqIO.parse(path, "fasta")]
    assert len(ids) == len(set(ids)), f"{protein}: duplicate IDs in pdb.fasta: {[i for i in ids if ids.count(i) > 1]}"


@pytest.mark.parametrize("protein", config.proteins())
def test_no_duplicate_alignment_ids(protein, data_root):
    path = data_root / protein / "refs" / "alignment.fasta"
    ids = [r.id for r in SeqIO.parse(path, "fasta")]
    assert len(ids) == len(set(ids)), f"{protein}: duplicate IDs in alignment.fasta"


# ---------------------------------------------------------------------------
# 4. Cross-file ID coverage
# ---------------------------------------------------------------------------

def test_query_and_pdb_have_same_ids(protein_fastas):
    protein, query, pdb, _ = protein_fastas
    query_ids = _query_ids(query)
    pdb_bare  = _pdb_ids(pdb)
    only_query = query_ids - pdb_bare
    only_pdb   = pdb_bare - query_ids
    assert not only_query, f"{protein}: IDs in query.fasta but not pdb.fasta: {only_query}"
    assert not only_pdb,   f"{protein}: IDs in pdb.fasta but not query.fasta: {only_pdb}"


def test_all_query_ids_in_alignment(protein_fastas):
    protein, query, _, align = protein_fastas
    missing = _query_ids(query) - set(align.keys())
    assert not missing, f"{protein}: query IDs absent from alignment.fasta: {missing}"


def test_all_pdb_ids_in_alignment(protein_fastas):
    protein, _, pdb, align = protein_fastas
    missing = set(pdb.keys()) - set(align.keys())
    assert not missing, f"{protein}: pdb IDs absent from alignment.fasta: {missing}"


def test_no_extra_ids_in_alignment(protein_fastas):
    protein, query, pdb, align = protein_fastas
    expected = _query_ids(query) | set(pdb.keys())
    extra = set(align.keys()) - expected
    assert not extra, f"{protein}: unexpected IDs in alignment.fasta: {extra}"


# ---------------------------------------------------------------------------
# 5. Alignment pair integrity
# ---------------------------------------------------------------------------

def test_alignment_pairs_same_length(protein_fastas):
    protein, query, _, align = protein_fastas
    mismatches = {}
    for id_ in _query_ids(query):
        pdb_key = id_ + "_pdb"
        if id_ not in align or pdb_key not in align:
            continue  # covered by coverage tests above
        q_len = len(align[id_])
        p_len = len(align[pdb_key])
        if q_len != p_len:
            mismatches[id_] = (q_len, p_len)
    assert not mismatches, f"{protein}: alignment row length mismatch: {mismatches}"


def test_alignment_rows_not_all_gaps(protein_fastas):
    protein, query, _, align = protein_fastas
    all_gaps = []
    for id_ in _query_ids(query):
        for key in (id_, id_ + "_pdb"):
            if key in align and set(align[key]) == {GAP}:
                all_gaps.append(key)
    assert not all_gaps, f"{protein}: alignment rows that are entirely gaps: {all_gaps}"


def test_alignment_pairs_have_common_columns(protein_fastas):
    protein, query, _, align = protein_fastas
    no_overlap = []
    for id_ in _query_ids(query):
        pdb_key = id_ + "_pdb"
        if id_ not in align or pdb_key not in align:
            continue
        q_seq, p_seq = align[id_], align[pdb_key]
        if len(q_seq) != len(p_seq):
            continue  # already caught above
        shared = any(q != GAP and p != GAP for q, p in zip(q_seq, p_seq))
        if not shared:
            no_overlap.append(id_)
    assert not no_overlap, f"{protein}: pairs with no overlapping non-gap column: {no_overlap}"
