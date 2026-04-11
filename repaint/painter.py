# repainting logic

import numpy as np


def build_new_types(labels, lengths, n_chains):
    """
    Build the flat new atom_type array for all chains.

    labels  : list of letters, e.g. ['A', 'B', 'A']
    lengths : list of ints,   e.g. [10, 60, 10]
    n_chains: int

    Letters are mapped to integers in order of first appearance (A→1, B→2, …).
    Returns a 1-D numpy int array of length n_chains * sum(lengths).
    """
    unique_labels = sorted(set(labels))
    letter_to_type = {label: i + 1 for i, label in enumerate(unique_labels)}

    one_chain = np.concatenate([
        np.full(length, letter_to_type[label], dtype=int)
        for label, length in zip(labels, lengths)
    ])

    return np.tile(one_chain, n_chains)


def repaint(df, labels, lengths):
    """
    Return a copy of df with atom_type replaced by the requested pattern.
    """
    n_chains = df['mol_id'].nunique()
    df = df.copy()
    df['atom_type'] = build_new_types(labels, lengths, n_chains)
    return df


def write_output(df, lines, atoms_start, velocities_start, output_path):
    """
    Write the repainted datafile, preserving the original header and
    everything from Velocities onward.
    """
    header_lines = lines[:atoms_start + 1]  # +1 keeps the blank line after "Atoms  # molecular"
    tail = lines[velocities_start - 1:]

    with open(output_path, 'w') as f:
        f.writelines(header_lines)
        df[['atom_id', 'mol_id', 'atom_type', 'x', 'y', 'z', 'ix', 'iy', 'iz']].to_csv(
            f, sep=' ', index=False, header=False
        )
        f.writelines(tail)
