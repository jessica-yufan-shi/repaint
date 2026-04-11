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
