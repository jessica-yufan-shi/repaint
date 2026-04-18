import textwrap
import numpy as np
import pytest
from repaint.painter import build_new_types, build_random_types, repaint, repaint_random, write_output
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
# build_random_types
# ---------------------------------------------------------------------------

def test_build_random_types_length():
    result = build_random_types(['A', 'B'], [0.4, 0.6], chain_length=80, n_chains=100, seed=0)
    assert len(result) == 80 * 100

def test_build_random_types_dtype():
    result = build_random_types(['A', 'B'], [0.5, 0.5], chain_length=10, n_chains=5, seed=0)
    assert result.dtype == int

def test_build_random_types_only_valid_types():
    result = build_random_types(['A', 'B', 'C'], [0.3, 0.3, 0.4], chain_length=50, n_chains=20, seed=0)
    assert set(result).issubset({1, 2, 3})

def test_build_random_types_seed_reproducibility():
    r1 = build_random_types(['A', 'B'], [0.4, 0.6], chain_length=80, n_chains=50, seed=42)
    r2 = build_random_types(['A', 'B'], [0.4, 0.6], chain_length=80, n_chains=50, seed=42)
    np.testing.assert_array_equal(r1, r2)

def test_build_random_types_different_seeds_differ():
    r1 = build_random_types(['A', 'B'], [0.5, 0.5], chain_length=80, n_chains=50, seed=1)
    r2 = build_random_types(['A', 'B'], [0.5, 0.5], chain_length=80, n_chains=50, seed=2)
    assert not np.array_equal(r1, r2)

def test_build_random_types_chains_are_independent():
    # With enough monomers, two chains sampled from seed 0 should not be identical
    result = build_random_types(['A', 'B'], [0.5, 0.5], chain_length=100, n_chains=10, seed=0)
    chains = result.reshape(10, 100)
    # Not all chains should be identical
    assert not all(np.array_equal(chains[0], chains[i]) for i in range(1, 10))

def test_build_random_types_composition_approximately_correct():
    # With many monomers, the overall A fraction should be close to 0.4
    result = build_random_types(['A', 'B'], [0.4, 0.6], chain_length=200, n_chains=500, seed=0)
    a_fraction = np.sum(result == 1) / len(result)
    assert abs(a_fraction - 0.4) < 0.02

def test_build_random_types_letter_mapping():
    # A→1, B→2 (alphabetical)
    result = build_random_types(['A', 'B'], [0.5, 0.5], chain_length=1000, n_chains=1, seed=0)
    assert set(result).issubset({1, 2})

def test_build_random_types_three_components_composition():
    result = build_random_types(['A', 'B', 'C'], [0.2, 0.5, 0.3], chain_length=200, n_chains=500, seed=7)
    total = len(result)
    assert abs(np.sum(result == 1) / total - 0.2) < 0.02
    assert abs(np.sum(result == 2) / total - 0.5) < 0.02
    assert abs(np.sum(result == 3) / total - 0.3) < 0.02


# ---------------------------------------------------------------------------
# repaint_random
# ---------------------------------------------------------------------------

def test_repaint_random_returns_copy(tmp_path):
    filepath = make_dummy_datafile(tmp_path, [1, 1, 2, 2], n_chains=10)
    df, *_ = load_datafile(filepath)
    original_types = df['atom_type'].tolist()
    repaint_random(df, ['A', 'B'], [0.5, 0.5], seed=0)
    assert df['atom_type'].tolist() == original_types

def test_repaint_random_other_columns_unchanged(tmp_path):
    filepath = make_dummy_datafile(tmp_path, [1, 1, 2, 2], n_chains=10)
    df, *_ = load_datafile(filepath)
    result = repaint_random(df, ['A', 'B'], [0.5, 0.5], seed=0)
    for col in ['atom_id', 'mol_id', 'x', 'y', 'z', 'ix', 'iy', 'iz']:
        assert result[col].tolist() == df[col].tolist()

def test_repaint_random_seed_reproducibility(tmp_path):
    filepath = make_dummy_datafile(tmp_path, [1, 1, 2, 2], n_chains=10)
    df, *_ = load_datafile(filepath)
    r1 = repaint_random(df, ['A', 'B'], [0.4, 0.6], seed=99)
    r2 = repaint_random(df, ['A', 'B'], [0.4, 0.6], seed=99)
    assert r1['atom_type'].tolist() == r2['atom_type'].tolist()


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
    # mapping is alphabetical (A→1, B→2) regardless of order in labels
    result = build_new_types(['B', 'A'], [2, 2], n_chains=1)
    expected = np.array([2, 2, 1, 1])
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
