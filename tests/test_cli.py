import textwrap
import pytest
from unittest.mock import patch
from repaint.cli import main


# ---------------------------------------------------------------------------
# Dummy datafile helper
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
# Argument validation (no file I/O needed)
# ---------------------------------------------------------------------------

def test_file_not_found(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--file", "nonexistent_file", "--pattern", "A-B", "--lengths", "5-5"])
    assert exc.value.code == 1
    assert "not found" in capsys.readouterr().out

def test_pattern_lengths_count_mismatch(tmp_path, capsys):
    filepath = make_dummy_datafile(tmp_path, [1, 1, 2, 2], n_chains=2)
    with pytest.raises(SystemExit) as exc:
        main(["--file", filepath, "--pattern", "A-B-A", "--lengths", "2-2"])
    assert exc.value.code == 1
    out = capsys.readouterr().out
    assert "pattern has 3 blocks" in out
    assert "lengths has 2 values" in out

def test_lengths_sum_mismatch(tmp_path, capsys):
    # chain length = 4, but lengths sum to 6
    filepath = make_dummy_datafile(tmp_path, [1, 1, 2, 2], n_chains=2)
    with pytest.raises(SystemExit) as exc:
        main(["--file", filepath, "--pattern", "A-B", "--lengths", "3-3"])
    assert exc.value.code == 1
    out = capsys.readouterr().out
    assert "lengths sum to 6" in out
    assert "chain length is 4" in out

def test_invalid_lengths_not_integers(tmp_path, capsys):
    filepath = make_dummy_datafile(tmp_path, [1, 2], n_chains=1)
    with pytest.raises(SystemExit) as exc:
        main(["--file", filepath, "--pattern", "A-B", "--lengths", "1-x"])
    assert exc.value.code == 1
    assert "integers" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# User says no
# ---------------------------------------------------------------------------

def test_user_aborts(tmp_path, capsys):
    filepath = make_dummy_datafile(tmp_path, [1, 1, 2, 2], n_chains=2)
    with patch("builtins.input", return_value="no"):
        with pytest.raises(SystemExit) as exc:
            main(["--file", filepath, "--pattern", "A-B", "--lengths", "2-2"])
    assert exc.value.code == 0
    assert "Aborted" in capsys.readouterr().out
    # no output file written
    assert not (tmp_path / "repainted_dummy_data").exists()


# ---------------------------------------------------------------------------
# Successful repaint
# ---------------------------------------------------------------------------

def test_successful_repaint_default_output(tmp_path, capsys):
    filepath = make_dummy_datafile(tmp_path, [1, 1, 2, 2], n_chains=2)
    with patch("builtins.input", return_value="yes"):
        with patch("repaint.cli.write_output") as mock_write:
            main(["--file", filepath, "--pattern", "A-B", "--lengths", "2-2"])
    mock_write.assert_called_once()
    out = capsys.readouterr().out
    assert "Written to: repainted_dummy_data" in out

def test_successful_repaint_custom_output(tmp_path, capsys):
    filepath = make_dummy_datafile(tmp_path, [1, 1, 2, 2], n_chains=2)
    out_path = str(tmp_path / "my_output")
    with patch("builtins.input", return_value="yes"):
        with patch("repaint.cli.write_output") as mock_write:
            main(["--file", filepath, "--pattern", "A-B", "--lengths", "2-2",
                  "--output", out_path])
    _, _, _, _, called_path = mock_write.call_args[0]
    assert called_path == out_path

def test_end_to_end_file_written(tmp_path, capsys):
    # Full pipeline without mocking write_output — verifies the file on disk
    filepath = make_dummy_datafile(tmp_path, [1, 1, 2, 2, 1, 1], n_chains=3)
    out_path = str(tmp_path / "out")
    with patch("builtins.input", return_value="yes"):
        main(["--file", filepath, "--pattern", "A-B-A", "--lengths", "2-2-2",
              "--output", out_path])

    from repaint.reader import load_datafile, get_chain_info
    import numpy as np
    df, *_ = load_datafile(out_path)
    n_chains, chain_length = get_chain_info(df)
    assert n_chains == 3
    assert chain_length == 6
    for mol, grp in df.groupby("mol_id"):
        np.testing.assert_array_equal(grp["atom_type"].values, [1, 1, 2, 2, 1, 1])
