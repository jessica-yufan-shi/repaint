import re
import numpy as np
import pandas as pd
from unittest.mock import patch
from colorama import Fore, Style
from repaint.display import build_color_map, render_chain, render_random_chain, show_confirmation, show_random_confirmation, save_composition_histogram

ansi_escape = re.compile(r'\x1b\[[0-9;]*m')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_repainted_df(labels, probs, chain_length, n_chains, seed=0):
    """Build a DataFrame as repaint_random would return it."""
    from repaint.painter import build_random_types
    types = build_random_types(labels, probs, chain_length, n_chains, seed=seed)
    mol_ids = np.repeat(np.arange(1, n_chains + 1), chain_length)
    atom_ids = np.arange(1, n_chains * chain_length + 1)
    return pd.DataFrame({
        'atom_id': atom_ids, 'mol_id': mol_ids, 'atom_type': types,
        'x': 0.0, 'y': 0.0, 'z': 0.0, 'ix': 0, 'iy': 0, 'iz': 0,
    })


# ---------------------------------------------------------------------------
# save_composition_histogram
# ---------------------------------------------------------------------------

def test_histogram_creates_png(tmp_path):
    df = make_repainted_df(['A', 'B'], [0.4, 0.6], chain_length=40, n_chains=50)
    out = str(tmp_path / "output")
    hist_path = save_composition_histogram(df, ['A', 'B'], [0.4, 0.6], out)
    assert hist_path == out + '_composition.png'
    assert (tmp_path / "output_composition.png").exists()

def test_histogram_three_components(tmp_path):
    df = make_repainted_df(['A', 'B', 'C'], [0.3, 0.4, 0.3], chain_length=60, n_chains=40)
    out = str(tmp_path / "output")
    hist_path = save_composition_histogram(df, ['A', 'B', 'C'], [0.3, 0.4, 0.3], out)
    assert (tmp_path / "output_composition.png").exists()

def test_histogram_returns_correct_path(tmp_path):
    df = make_repainted_df(['A', 'B'], [0.5, 0.5], chain_length=20, n_chains=10)
    out = str(tmp_path / "my_output")
    result = save_composition_histogram(df, ['A', 'B'], [0.5, 0.5], out)
    assert result == str(tmp_path / "my_output_composition.png")

def test_histogram_single_component(tmp_path):
    # Edge case: single component (100% A)
    df = make_repainted_df(['A'], [1.0], chain_length=20, n_chains=10)
    out = str(tmp_path / "output")
    hist_path = save_composition_histogram(df, ['A'], [1.0], out)
    assert (tmp_path / "output_composition.png").exists()


# ---------------------------------------------------------------------------
# build_color_map
# ---------------------------------------------------------------------------

def test_color_map_consistent_across_patterns():
    # A and B appear in both; color for A must be the same in both
    color_map = build_color_map(['A', 'B', 'A'], ['A', 'B', 'A'])
    assert color_map['A'] == color_map['A']  # trivially true, but explicit
    assert color_map['A'] != color_map['B']

def test_color_map_order_of_appearance():
    # C only appears in the second list; it should still get a slot after A and B
    color_map = build_color_map(['A', 'B'], ['A', 'B', 'C'])
    from repaint.display import _PALETTE
    assert color_map['A'] == _PALETTE[0]
    assert color_map['B'] == _PALETTE[1]
    assert color_map['C'] == _PALETTE[2]

def test_color_map_no_duplicates():
    color_map = build_color_map(['A', 'B', 'A', 'B'], ['A'])
    assert len(color_map) == 2   # only A and B


# ---------------------------------------------------------------------------
# render_chain
# ---------------------------------------------------------------------------

def test_render_chain_contains_letters():
    color_map = build_color_map(['A', 'B'], ['A', 'B'])
    output = render_chain(['A', 'B', 'A'], [2, 4, 2], color_map)
    assert 'AA' in output
    assert 'BBBB' in output

def test_render_chain_brackets():
    color_map = build_color_map(['A'], ['A'])
    output = render_chain(['A'], [3], color_map)
    # strip ANSI and check brackets
    plain = output.replace(Style.RESET_ALL, '').replace(color_map['A'], '')
    assert plain == '[AAA]'

def test_render_chain_aba_structure():
    import re
    ansi_escape = re.compile(r'\x1b\[[0-9;]*m')
    color_map = build_color_map(['A', 'B'], ['A', 'B'])
    output = render_chain(['A', 'B', 'A'], [1, 2, 1], color_map)
    plain = ansi_escape.sub('', output)
    assert plain == '[A][BB][A]'

def test_render_chain_color_applied():
    color_map = build_color_map(['A', 'B'], ['A', 'B'])
    output = render_chain(['A', 'B'], [2, 3], color_map)
    assert color_map['A'] in output
    assert color_map['B'] in output


# ---------------------------------------------------------------------------
# show_confirmation
# ---------------------------------------------------------------------------

def test_show_confirmation_yes(capsys):
    with patch('builtins.input', return_value='yes'):
        result = show_confirmation(
            filename='test_data',
            n_chains=10,
            chain_length=20,
            current_labels=['A', 'B', 'A'],
            current_lengths=[5, 10, 5],
            requested_labels=['A', 'B', 'A'],
            requested_lengths=[4, 12, 4],
        )
    assert result is True
    out = capsys.readouterr().out
    assert 'File: test_data' in out
    assert '10 chains' in out
    assert '20 monomers' in out

def test_show_confirmation_no(capsys):
    with patch('builtins.input', return_value='no'):
        result = show_confirmation(
            filename='test_data',
            n_chains=5,
            chain_length=10,
            current_labels=['A', 'B'],
            current_lengths=[5, 5],
            requested_labels=['A', 'B', 'A'],
            requested_lengths=[2, 6, 2],
        )
    assert result is False

# ---------------------------------------------------------------------------
# render_random_chain
# ---------------------------------------------------------------------------

def test_render_random_chain_contains_only_valid_letters():
    color_map = build_color_map(['A', 'B'], ['A', 'B'])
    type_to_label = {1: 'A', 2: 'B'}
    types = np.array([1, 2, 1, 1, 2])
    output = render_random_chain(types, type_to_label, color_map)
    plain = ansi_escape.sub('', output)
    assert set(plain) == {'A', 'B'}

def test_render_random_chain_correct_length():
    color_map = build_color_map(['A', 'B'], ['A', 'B'])
    type_to_label = {1: 'A', 2: 'B'}
    types = np.array([1, 2, 1, 2, 1])
    output = render_random_chain(types, type_to_label, color_map)
    plain = ansi_escape.sub('', output)
    assert len(plain) == 5

def test_render_random_chain_truncates():
    color_map = build_color_map(['A', 'B'], ['A', 'B'])
    type_to_label = {1: 'A', 2: 'B'}
    types = np.ones(100, dtype=int)
    output = render_random_chain(types, type_to_label, color_map, max_width=10)
    plain = ansi_escape.sub('', output)
    assert plain.endswith('...')
    assert len(plain) == 13  # 10 letters + '...'

def test_render_random_chain_no_truncation_when_short():
    color_map = build_color_map(['A', 'B'], ['A', 'B'])
    type_to_label = {1: 'A', 2: 'B'}
    types = np.array([1, 2, 1])
    output = render_random_chain(types, type_to_label, color_map, max_width=10)
    plain = ansi_escape.sub('', output)
    assert '...' not in plain


# ---------------------------------------------------------------------------
# show_random_confirmation
# ---------------------------------------------------------------------------

def test_show_random_confirmation_yes(capsys):
    with patch('builtins.input', return_value='yes'):
        result = show_random_confirmation(
            filename='test_data', n_chains=50, chain_length=80,
            current_labels=['A', 'B'], current_lengths=[40, 40],
            labels=['A', 'B'], probs=[0.4, 0.6], seed=0,
        )
    assert result is True

def test_show_random_confirmation_no(capsys):
    with patch('builtins.input', return_value='no'):
        result = show_random_confirmation(
            filename='test_data', n_chains=50, chain_length=80,
            current_labels=['A', 'B'], current_lengths=[40, 40],
            labels=['A', 'B'], probs=[0.4, 0.6], seed=0,
        )
    assert result is False

def test_show_random_confirmation_output_structure(capsys):
    with patch('builtins.input', return_value='no'):
        show_random_confirmation(
            filename='myfile', n_chains=100, chain_length=80,
            current_labels=['A', 'B'], current_lengths=[40, 40],
            labels=['A', 'B'], probs=[0.4, 0.6], seed=42,
        )
    out = capsys.readouterr().out
    assert 'File: myfile' in out
    assert '100 chains' in out
    assert '80 monomers' in out
    assert 'random' in out
    assert 'A: 40.0%' in out
    assert 'B: 60.0%' in out
    assert 'Sample chains (3 of 100)' in out

def test_show_random_confirmation_sample_chains_reproducible(capsys):
    kwargs = dict(
        filename='f', n_chains=50, chain_length=40,
        current_labels=['A', 'B'], current_lengths=[20, 20],
        labels=['A', 'B'], probs=[0.5, 0.5], seed=7,
    )
    with patch('builtins.input', return_value='no'):
        show_random_confirmation(**kwargs)
    out1 = capsys.readouterr().out

    with patch('builtins.input', return_value='no'):
        show_random_confirmation(**kwargs)
    out2 = capsys.readouterr().out

    assert out1 == out2


# ---------------------------------------------------------------------------
# show_confirmation
# ---------------------------------------------------------------------------

def test_show_confirmation_pattern_lines(capsys):
    with patch('builtins.input', return_value='no'):
        show_confirmation(
            filename='myfile',
            n_chains=100,
            chain_length=80,
            current_labels=['A', 'B', 'A', 'C'],
            current_lengths=[5, 60, 5, 10],
            requested_labels=['A', 'B', 'A'],
            requested_lengths=[10, 60, 10],
        )
    out = capsys.readouterr().out
    assert 'Current pattern:  A-B-A-C  |  lengths: 5-60-5-10' in out
    assert 'Requested pattern:  A-B-A  |  lengths: 10-60-10' in out
