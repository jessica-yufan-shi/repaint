import textwrap
import numpy as np
import pytest
from repaint.painter import build_new_types, repaint, write_output
from repaint.reader import load_datafile


# ---------------------------------------------------------------------------
# Helpers (shared with test_reader — duplicated here to keep tests isolated)
# ---------------------------------------------------------------------------

def make_dummy_datafile(tmp_path, pattern_types, n_chains=2):
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
# build_new_types
# ---------------------------------------------------------------------------

def test_build_new_types_length():
    result = build_new_types(['A', 'B', 'A'], [2, 4, 2], n_chains=3)
    assert len(result) == (2 + 4 + 2) * 3

def test_build_new_types_single_chain_values():
    result = build_new_types(['A', 'B', 'A'], [2, 3, 2], n_chains=1)
    expected = np.array([1, 1, 2, 2, 2, 1, 1])
    np.testing.assert_array_equal(result, expected)

def test_build_new_types_tiled():
    result = build_new_types(['A', 'B'], [2, 2], n_chains=3)
    expected = np.array([1, 1, 2, 2,  1, 1, 2, 2,  1, 1, 2, 2])
    np.testing.assert_array_equal(result, expected)

def test_build_new_types_letter_mapping():
    # B appears first in labels → B gets type 1, A gets type 2
    result = build_new_types(['B', 'A'], [2, 2], n_chains=1)
    expected = np.array([1, 1, 2, 2])
    np.testing.assert_array_equal(result, expected)

def test_build_new_types_repeated_letter():
    # A-B-A: A maps to 1, B maps to 2
    result = build_new_types(['A', 'B', 'A'], [1, 2, 1], n_chains=1)
    expected = np.array([1, 2, 2, 1])
    np.testing.assert_array_equal(result, expected)

def test_build_new_types_four_block():
    result = build_new_types(['A', 'B', 'C', 'D'], [1, 1, 1, 1], n_chains=1)
    np.testing.assert_array_equal(result, np.array([1, 2, 3, 4]))


# ---------------------------------------------------------------------------
# repaint
# ---------------------------------------------------------------------------

def test_repaint_returns_copy(tmp_path):
    filepath = make_dummy_datafile(tmp_path, [1, 1, 2, 2], n_chains=2)
    df, *_ = load_datafile(filepath)
    original_types = df['atom_type'].tolist()
    repainted = repaint(df, ['A', 'B'], [2, 2])
    # original df unchanged
    assert df['atom_type'].tolist() == original_types

def test_repaint_correct_types(tmp_path):
    filepath = make_dummy_datafile(tmp_path, [1, 1, 2, 2], n_chains=2)
    df, *_ = load_datafile(filepath)
    repainted = repaint(df, ['A', 'B', 'A'], [1, 2, 1])
    # each chain should be [1, 2, 2, 1]
    for mol, grp in repainted.groupby('mol_id'):
        np.testing.assert_array_equal(
            grp['atom_type'].values, [1, 2, 2, 1]
        )

def test_repaint_other_columns_unchanged(tmp_path):
    filepath = make_dummy_datafile(tmp_path, [1, 1, 2, 2, 1], n_chains=2)
    df, *_ = load_datafile(filepath)
    repainted = repaint(df, ['A', 'B', 'A'], [2, 1, 2])
    for col in ['atom_id', 'mol_id', 'x', 'y', 'z', 'ix', 'iy', 'iz']:
        assert repainted[col].tolist() == df[col].tolist()


# ---------------------------------------------------------------------------
# write_output
# ---------------------------------------------------------------------------

def test_write_output_roundtrip(tmp_path):
    pattern = [1, 1, 2, 2, 1, 1]
    filepath = make_dummy_datafile(tmp_path, pattern, n_chains=2)
    df, lines, atoms_start, velocities_start = load_datafile(filepath)

    repainted = repaint(df, ['A', 'B', 'A'], [2, 2, 2])
    out_path = str(tmp_path / "out_data")
    write_output(repainted, lines, atoms_start, velocities_start, out_path)

    # Read back and verify atom types
    df2, *_ = load_datafile(out_path)
    for mol, grp in df2.groupby('mol_id'):
        np.testing.assert_array_equal(grp['atom_type'].values, [1, 1, 2, 2, 1, 1])

def test_write_output_preserves_velocities(tmp_path):
    pattern = [1, 2, 1]
    filepath = make_dummy_datafile(tmp_path, pattern, n_chains=2)
    df, lines, atoms_start, velocities_start = load_datafile(filepath)

    out_path = str(tmp_path / "out_data")
    write_output(df, lines, atoms_start, velocities_start, out_path)

    with open(out_path) as f:
        content = f.read()
    assert 'Velocities' in content

def test_write_output_preserves_header(tmp_path):
    pattern = [1, 2]
    filepath = make_dummy_datafile(tmp_path, pattern, n_chains=1)
    df, lines, atoms_start, velocities_start = load_datafile(filepath)

    out_path = str(tmp_path / "out_data")
    write_output(df, lines, atoms_start, velocities_start, out_path)

    with open(out_path) as f:
        first_line = f.readline()
    assert first_line == lines[0]
