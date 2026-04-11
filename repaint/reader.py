# file reading and chain detection

import pandas as pd


def find_section_boundaries(lines):
    atoms_start = next(i for i, l in enumerate(lines) if l.strip().replace('  ', ' ') == 'Atoms # molecular') + 1
    velocities_start = next(i for i, l in enumerate(lines) if l.strip() == 'Velocities')
    return atoms_start, velocities_start


def load_datafile(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    atoms_start, velocities_start = find_section_boundaries(lines)

    headers = ['atom_id', 'mol_id', 'atom_type', 'x', 'y', 'z', 'ix', 'iy', 'iz']
    df = pd.read_csv(
        filepath,
        names=headers,
        skiprows=atoms_start,
        nrows=velocities_start - atoms_start - 2,
        delimiter=' ',
    )

    # Step 2b: sort by mol_id then atom_id
    df = df.sort_values(['mol_id', 'atom_id']).reset_index(drop=True)

    return df, lines, atoms_start, velocities_start


def check_contiguity(df):
    def is_contiguous(ids):
        return ids.max() - ids.min() + 1 == len(ids)

    bad = df.groupby('mol_id')['atom_id'].apply(is_contiguous)
    if not bad.all():
        print(
            "Warning: atom IDs within some chains are not contiguous after grouping "
            "by mol_id. File may be corrupt or have unexpected format."
        )


def get_chain_info(df):
    n_chains = df['mol_id'].nunique()
    chain_length = len(df) // n_chains
    return n_chains, chain_length


def _rle(sequence):
    """Run-length encode a list. Returns list of (value, count) tuples."""
    result = []
    if not sequence:
        return result
    current, count = sequence[0], 1
    for val in sequence[1:]:
        if val == current:
            count += 1
        else:
            result.append((current, count))
            current, count = val, 1
    result.append((current, count))
    return result


def detect_current_pattern(df):
    """
    Infer block pattern from the first chain via run-length encoding.

    Returns:
        labels       : list of letters, e.g. ['A', 'B', 'A']
        lengths      : list of ints,   e.g. [3, 4, 3]
        type_to_letter: dict mapping int atom type -> letter
    """
    first_mol = df['mol_id'].iloc[0]
    first_chain = df[df['mol_id'] == first_mol].sort_values('atom_id')
    types = first_chain['atom_type'].tolist()

    encoded = _rle(types)

    type_to_letter = {}
    next_code = ord('A')
    labels, lengths = [], []
    for typ, count in encoded:
        if typ not in type_to_letter:
            type_to_letter[typ] = chr(next_code)
            next_code += 1
        labels.append(type_to_letter[typ])
        lengths.append(count)

    return labels, lengths, type_to_letter
