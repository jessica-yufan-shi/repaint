# argument parsing and main entry point

import argparse
import os
import sys

from colorama import Fore, Style, init
from repaint.reader import load_datafile, check_contiguity, get_chain_info, detect_current_pattern
from repaint.display import show_confirmation
from repaint.painter import repaint, write_output

init(autoreset=True)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Repaint copolymer atom types in a LAMMPS datafile."
    )
    parser.add_argument("--file", required=True, help="Path to the LAMMPS datafile")
    parser.add_argument("--pattern", required=True, help="Block pattern, e.g. A-B-A")
    parser.add_argument("--lengths", required=True, help="Block lengths, e.g. 10-60-10")
    parser.add_argument(
        "--output", default=None,
        help="Output filename (default: repainted_<input filename>)"
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # Validate file exists
    if not os.path.isfile(args.file):
        print(f"Error: file '{args.file}' not found.")
        sys.exit(1)

    # Parse pattern and lengths
    labels = args.pattern.split("-")
    try:
        lengths = [int(x) for x in args.lengths.split("-")]
    except ValueError:
        print("Error: lengths must be dash-separated integers, e.g. 10-60-10.")
        sys.exit(1)

    # Validate count match before reading the (potentially large) file
    if len(labels) != len(lengths):
        print(
            f"Error: pattern has {len(labels)} blocks but lengths has {len(lengths)} values."
        )
        sys.exit(1)

    # Read file
    df, lines, atoms_start, velocities_start = load_datafile(args.file)
    check_contiguity(df)
    n_chains, chain_length = get_chain_info(df)

    # Validate lengths sum matches chain length
    if sum(lengths) != chain_length:
        print(
            f"Error: lengths sum to {sum(lengths)} but chain length is {chain_length}."
        )
        sys.exit(1)

    # Detect current pattern from file
    current_labels, current_lengths, _ = detect_current_pattern(df)

    # Determine output path
    os.makedirs("results", exist_ok=True)
    filename = os.path.basename(args.file)
    if args.output:
        output_path = args.output
    else:
        pattern_tag = f"{args.pattern}_{args.lengths}"
        output_path = os.path.join("results", f"repainted_{filename}_{pattern_tag}")

    # Warn if requested pattern introduces new atom types
    original_n_types = df['atom_type'].nunique()
    requested_n_types = len(set(labels))
    if requested_n_types > original_n_types:
        print(
            f"{Fore.YELLOW}Warning: atom types increased from {original_n_types} to {requested_n_types}. "
            f"Remember to update Atom Types, Masses, and Pair Coeffs (and any other forcefield "
            f"sections) in the output file before running LAMMPS.{Style.RESET_ALL}"
        )

    # Show confirmation preview
    confirmed = show_confirmation(
        filename=filename,
        n_chains=n_chains,
        chain_length=chain_length,
        current_labels=current_labels,
        current_lengths=current_lengths,
        requested_labels=labels,
        requested_lengths=lengths,
    )

    if not confirmed:
        print("Aborted. No file written.")
        sys.exit(0)

    # Repaint and write
    df_repainted = repaint(df, labels, lengths)
    write_output(df_repainted, lines, atoms_start, velocities_start, output_path)
    print(f"Written to: {output_path}")
