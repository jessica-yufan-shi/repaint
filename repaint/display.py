# confirmation preview and colored output

from colorama import Fore, Style, init

init(autoreset=True)

_PALETTE = [
    Fore.CYAN,
    Fore.GREEN,
    Fore.YELLOW,
    Fore.MAGENTA,
    Fore.RED,
    Fore.BLUE,
    Fore.WHITE,
]


def build_color_map(labels_a, labels_b):
    """
    Assign a color to each unique letter seen across both label lists.
    Order of first appearance determines the color slot, so A is always
    the same color regardless of which pattern it appears in.
    """
    seen = []
    for label in labels_a + labels_b:
        if label not in seen:
            seen.append(label)
    return {letter: _PALETTE[i % len(_PALETTE)] for i, letter in enumerate(seen)}


def render_chain(labels, lengths, color_map):
    """Return a colored string like [AAAA][BBBBBB][AAAA]."""
    parts = []
    for label, length in zip(labels, lengths):
        color = color_map.get(label, Fore.WHITE)
        parts.append(f"{color}[{label * length}]{Style.RESET_ALL}")
    return "".join(parts)


def render_random_chain(types, type_to_label, color_map, max_width=80):
    """
    Render a random chain as a colored sequence of individual letters.
    Truncates to max_width visible characters and appends '...' if longer.
    """
    display_types = types[:max_width]
    parts = []
    for t in display_types:
        label = type_to_label[t]
        color = color_map.get(label, Fore.WHITE)
        parts.append(f"{color}{label}{Style.RESET_ALL}")
    result = "".join(parts)
    if len(types) > max_width:
        result += "..."
    return result


def show_random_confirmation(
    filename,
    n_chains,
    chain_length,
    current_labels,
    current_lengths,
    labels,
    probs,
    seed,
):
    """
    Print the random-mode confirmation preview and prompt the user.
    Shows the current pattern, target composition, and 3 sample chains.
    Returns True if the user types 'yes', False otherwise.
    """
    from repaint.painter import build_random_types

    color_map = build_color_map(current_labels, labels)
    type_to_label = {i + 1: label for i, label in enumerate(sorted(set(labels)))}

    current_pattern_str = "-".join(current_labels)
    current_lengths_str = "-".join(str(l) for l in current_lengths)

    print(f"File: {filename}")
    print(f"Chains: {n_chains} chains, {chain_length} monomers each")
    print()
    print(f"Current pattern:  {current_pattern_str}  |  lengths: {current_lengths_str}")
    print(render_chain(current_labels, current_lengths, color_map))
    print()

    comp_str = "  ".join(f"{label}: {prob * 100:.1f}%" for label, prob in zip(labels, probs))
    print(f"Requested: random  |  composition: {comp_str}")

    preview = build_random_types(labels, probs, chain_length, n_chains=3, seed=seed)
    preview = preview.reshape(3, chain_length)

    print(f"Sample chains (3 of {n_chains}):")
    for chain_types in preview:
        print(f"  {render_random_chain(chain_types, type_to_label, color_map)}")
    print()

    answer = input("Repaint? [yes/no]: ").strip().lower()
    return answer == "yes"


_PLOT_COLORS = [
    '#8D9EE3', '#E16F6E', '#C6D034',
    'mediumpurple', 'sandybrown', 'lightcoral', 'turquoise', 'khaki',
]


def save_composition_histogram(df_repainted, labels, probs, output_path):
    """
    Compute per-chain fractions for each component and save a histogram as a PNG.

    One subplot per component. Each subplot shows the distribution of per-chain
    fractions across all chains, a dashed line at the target, and a dotted line
    at the actual mean. Saved as <output_path>_composition.png.

    Returns the path to the saved PNG.
    """
    import math
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    unique_labels = sorted(set(labels))
    letter_to_type = {label: i + 1 for i, label in enumerate(unique_labels)}

    # Per-chain type counts → fractions
    chain_counts = df_repainted.groupby(['mol_id', 'atom_type']).size().unstack(fill_value=0)
    chain_lengths = chain_counts.sum(axis=1)
    chain_fractions = chain_counts.divide(chain_lengths, axis=0)

    n_components = len(labels)
    n_chains = len(chain_lengths)
    chain_length = int(chain_lengths.iloc[0])
    n_bins = min(30, max(10, n_chains // 5))

    ncols = min(n_components, 4)
    nrows = math.ceil(n_components / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows), sharey=False)
    axes = [axes] if n_components == 1 else list(np.array(axes).ravel())

    label_to_target = dict(zip(labels, probs))

    for i, label in enumerate(labels):
        ax = axes[i]
        atom_type = letter_to_type[label]
        target = label_to_target[label]

        frac = chain_fractions[atom_type].values if atom_type in chain_fractions.columns else [0.0] * n_chains
        mean = sum(frac) / len(frac)
        variance = sum((f - mean) ** 2 for f in frac) / len(frac)
        std = variance ** 0.5

        color = _PLOT_COLORS[i % len(_PLOT_COLORS)]
        ax.hist(frac, bins=n_bins, color=color, edgecolor='none', alpha=0.85)
        ax.axvline(target, color='black',   linestyle='--', linewidth=1.5, label=f'target {target:.1%}')
        ax.axvline(mean,   color='dimgray', linestyle=':',  linewidth=1.5, label=f'mean {mean:.1%}')
        ax.set_xlim(0, 1)
        ax.set_xlabel(f'Fraction of {label}')
        ax.set_ylabel('Number of chains')
        ax.set_title(f'{label}  (target {target:.1%}, mean {mean:.1%}, σ {std:.1%})')
        ax.legend(fontsize=8)

    for ax in axes[n_components:]:
        ax.set_visible(False)

    fig.suptitle(
        f'Per-chain composition  ({n_chains} chains × {chain_length} monomers)',
        y=1.02,
    )
    fig.tight_layout()

    hist_path = output_path + '_composition.png'
    fig.savefig(hist_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

    return hist_path


def show_confirmation(
    filename,
    n_chains,
    chain_length,
    current_labels,
    current_lengths,
    requested_labels,
    requested_lengths,
):
    """
    Print the side-by-side preview and prompt the user.
    Returns True if the user types 'yes', False otherwise.
    """
    color_map = build_color_map(current_labels, requested_labels)

    current_pattern_str = "-".join(current_labels)
    current_lengths_str = "-".join(str(l) for l in current_lengths)
    requested_pattern_str = "-".join(requested_labels)
    requested_lengths_str = "-".join(str(l) for l in requested_lengths)

    print(f"File: {filename}")
    print(f"Chains: {n_chains} chains, {chain_length} monomers each")
    print()
    print(f"Current pattern:  {current_pattern_str}  |  lengths: {current_lengths_str}")
    print(render_chain(current_labels, current_lengths, color_map))
    print()
    print(f"Requested pattern:  {requested_pattern_str}  |  lengths: {requested_lengths_str}")
    print(render_chain(requested_labels, requested_lengths, color_map))
    print()

    answer = input("Repaint? [yes/no]: ").strip().lower()
    return answer == "yes"
