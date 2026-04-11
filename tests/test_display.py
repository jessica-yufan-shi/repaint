from unittest.mock import patch
from colorama import Fore, Style
from repaint.display import build_color_map, render_chain, show_confirmation


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
