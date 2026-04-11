import textwrap
import pytest
import pandas as pd
from repaint.reader import (
    find_section_boundaries,
    load_datafile,
    check_contiguity,
    get_chain_info,
    detect_current_pattern,
    _rle,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_dummy_datafile(tmp_path, pattern_types, n_chains=2):
    """
    Write a minimal LAMMPS datafile with the given per-chain atom type list.
    pattern_types: list of ints, e.g. [1,1,2,2,2,1,1] — one chain worth
    n_chains     : number of identical chains to tile
    """
    chain_len = len(pattern_types)
    n_atoms = chain_len * n_chains
    n_types = len(set(pattern_types))

    header = textwrap.dedent(f"""\
        # dummy LAMMPS datafile

        {n_atoms} atoms
        {n_types} atom types

        0.0 10.0 xlo xhi
        0.0 10.0 ylo yhi
        0.0 10.0 zlo zhi

        Masses

        """)
    for t in range(1, n_types + 1):
        header += f"{t} 1.0\n"
    header += "\nAtoms  # molecular\n\n"

    atom_lines = []
    atom_id = 1
    for mol in range(1, n_chains + 1):
        for pos, typ in enumerate(pattern_types):
            atom_lines.append(
                f"{atom_id} {mol} {typ} {float(pos)} 0.0 0.0 0 0 0\n"
            )
            atom_id += 1

    tail = "\nVelocities\n\n"
    for i in range(1, n_atoms + 1):
        tail += f"{i} 0.0 0.0 0.0\n"

    p = tmp_path / "dummy_data"
    p.write_text(header + "".join(atom_lines) + tail)
    return str(p)


# ---------------------------------------------------------------------------
# _rle
# ---------------------------------------------------------------------------

def test_rle_basic():
    assert _rle([1, 1, 2, 2, 2, 1]) == [(1, 2), (2, 3), (1, 1)]

def test_rle_single():
    assert _rle([3]) == [(3, 1)]

def test_rle_empty():
    assert _rle([]) == []

def test_rle_all_same():
    assert _rle([2, 2, 2]) == [(2, 3)]


# ---------------------------------------------------------------------------
# find_section_boundaries
# ---------------------------------------------------------------------------

def test_find_section_boundaries(tmp_path):
    pattern = [1, 1, 2, 2, 1, 1]
    filepath = make_dummy_datafile(tmp_path, pattern, n_chains=1)
    with open(filepath) as f:
        lines = f.readlines()
    atoms_start, velocities_start = find_section_boundaries(lines)
    assert lines[atoms_start - 1].strip() == 'Atoms  # molecular'
    assert lines[velocities_start].strip() == 'Velocities'


# ---------------------------------------------------------------------------
# load_datafile
# ---------------------------------------------------------------------------

def test_load_datafile_shape(tmp_path):
    pattern = [1, 1, 2, 2, 2, 1, 1]
    filepath = make_dummy_datafile(tmp_path, pattern, n_chains=3)
    df, lines, atoms_start, velocities_start = load_datafile(filepath)
    assert len(df) == len(pattern) * 3
    assert list(df.columns) == ['atom_id', 'mol_id', 'atom_type', 'x', 'y', 'z', 'ix', 'iy', 'iz']

def test_load_datafile_sorted(tmp_path):
    pattern = [1, 2, 1]
    filepath = make_dummy_datafile(tmp_path, pattern, n_chains=2)
    df, *_ = load_datafile(filepath)
    # must be sorted by mol_id then atom_id
    assert (df['mol_id'].diff().fillna(0) >= 0).all()
    for mol, grp in df.groupby('mol_id'):
        assert (grp['atom_id'].diff().fillna(1) > 0).all()


# ---------------------------------------------------------------------------
# get_chain_info
# ---------------------------------------------------------------------------

def test_get_chain_info(tmp_path):
    pattern = [1, 1, 2, 2, 2]
    filepath = make_dummy_datafile(tmp_path, pattern, n_chains=4)
    df, *_ = load_datafile(filepath)
    n_chains, chain_length = get_chain_info(df)
    assert n_chains == 4
    assert chain_length == 5


# ---------------------------------------------------------------------------
# check_contiguity
# ---------------------------------------------------------------------------

def test_check_contiguity_clean(tmp_path, capsys):
    pattern = [1, 2, 1]
    filepath = make_dummy_datafile(tmp_path, pattern, n_chains=2)
    df, *_ = load_datafile(filepath)
    check_contiguity(df)
    assert capsys.readouterr().out == ''

def test_check_contiguity_warns(capsys):
    # manually construct a df with a gap in atom_ids for mol 1
    df = pd.DataFrame({
        'atom_id':   [1, 3, 5, 4, 6, 7],   # mol 1 has a gap: 1,3,5
        'mol_id':    [1, 1, 1, 2, 2, 2],
        'atom_type': [1, 2, 1, 1, 2, 1],
        'x': [0]*6, 'y': [0]*6, 'z': [0]*6,
        'ix': [0]*6, 'iy': [0]*6, 'iz': [0]*6,
    })
    check_contiguity(df)
    assert 'Warning' in capsys.readouterr().out


# ---------------------------------------------------------------------------
# detect_current_pattern
# ---------------------------------------------------------------------------

def test_detect_pattern_aba(tmp_path):
    pattern = [1, 1, 2, 2, 2, 2, 1, 1]   # A(2)-B(4)-A(2)
    filepath = make_dummy_datafile(tmp_path, pattern, n_chains=2)
    df, *_ = load_datafile(filepath)
    labels, lengths, type_to_letter = detect_current_pattern(df)
    assert labels == ['A', 'B', 'A']
    assert lengths == [2, 4, 2]
    assert type_to_letter == {1: 'A', 2: 'B'}

def test_detect_pattern_abac(tmp_path):
    pattern = [1]*5 + [2]*10 + [3]*5 + [4]*5
    filepath = make_dummy_datafile(tmp_path, pattern, n_chains=1)
    df, *_ = load_datafile(filepath)
    labels, lengths, type_to_letter = detect_current_pattern(df)
    assert labels == ['A', 'B', 'C', 'D']
    assert lengths == [5, 10, 5, 5]
