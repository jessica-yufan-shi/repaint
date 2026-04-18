# argument parsing and main entry point

import argparse
import os
import sys

from colorama import Fore, Style, init
from repaint.reader import load_datafile, check_contiguity, get_chain_info, detect_current_pattern, parse_composition
from repaint.display import show_confirmation, show_random_confirmation, save_composition_histogram
from repaint.painter import repaint, repaint_random, write_output

init(autoreset=True)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Repaint copolymer atom types in a LAMMPS datafile."
    )
    parser.add_argument("--file", required=True, help="Path to the LAMMPS datafile")
    parser.add_argument(
        "--mode", choices=["block", "random"], default="block",
        help="Repainting mode: 'block' (default) or 'random'"
    )
    # block mode arguments
    parser.add_argument("--pattern", default=None, help="Block pattern, e.g. A-B-A (--mode block)")
    parser.add_argument("--lengths", default=None, help="Block lengths, e.g. 10-60-10 (--mode block)")
    # random mode arguments
    parser.add_argument("--composition", default=None, help="Composition, e.g. A:0.4-B:0.6 (--mode random)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility (--mode random)")
    parser.add_argument(
        "--output", default=None,
        help="Output filename (default: repainted_<input filename>)"
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # Validate mode-specific arguments
    if args.mode == "block":
        if not args.pattern:
            print("Error: --pattern is required for --mode block.")
            sys.exit(1)
        if not args.lengths:
            print("Error: --lengths is required for --mode block.")
            sys.exit(1)
        if args.composition is not None:
            print("Error: --composition is only valid for --mode random.")
            sys.exit(1)
        if args.seed is not None:
            print("Error: --seed is only valid for --mode random.")
            sys.exit(1)
    elif args.mode == "random":
        if args.composition is None:
            print("Error: --composition is required for --mode random.")
            sys.exit(1)
        if args.seed is None:
            print("Error: --seed is required for --mode random.")
            sys.exit(1)
        if args.pattern is not None:
            print("Error: --pattern is not valid for --mode random.")
            sys.exit(1)
        if args.lengths is not None:
            print("Error: --lengths is not valid for --mode random.")
            sys.exit(1)

    # Validate file exists
    if not os.path.isfile(args.file):
        print(f"Error: file '{args.file}' not found.")
        sys.exit(1)

    # Parse mode-specific inputs
    if args.mode == "block":
        labels = args.pattern.split("-")
        try:
            lengths = [int(x) for x in args.lengths.split("-")]
        except ValueError:
            print("Error: lengths must be dash-separated integers, e.g. 10-60-10.")
            sys.exit(1)

        if len(labels) != len(lengths):
            print(
                f"Error: pattern has {len(labels)} blocks but lengths has {len(lengths)} values."
            )
            sys.exit(1)

    else:  # random
        try:
            composition = parse_composition(args.composition)
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)
        labels = [label for label, _ in composition]
        probs  = [prob  for _, prob  in composition]

    # Read file (shared)
    df, lines, atoms_start, velocities_start = load_datafile(args.file)
    check_contiguity(df)
    n_chains, chain_length = get_chain_info(df)

    # Block-only: validate lengths sum
    if args.mode == "block" and sum(lengths) != chain_length:
        print(
            f"Error: lengths sum to {sum(lengths)} but chain length is {chain_length}."
        )
        sys.exit(1)

    # Detect current pattern (shared)
    current_labels, current_lengths, _ = detect_current_pattern(df)

    # Determine output path
    os.makedirs("results", exist_ok=True)
    filename = os.path.basename(args.file)
    if args.output:
        output_path = args.output
    elif args.mode == "block":
        pattern_tag = f"{args.pattern}_{args.lengths}"
        output_path = os.path.join("results", f"repainted_{filename}_{pattern_tag}")
    else:
        composition_tag = args.composition.replace(":", "")
        output_path = os.path.join("results", f"repainted_{filename}_random_{composition_tag}_seed{args.seed}")

    # Warn if requested pattern introduces new atom types (shared)
    original_n_types = df['atom_type'].nunique()
    requested_n_types = len(set(labels))
    if requested_n_types > original_n_types:
        print(
            f"{Fore.YELLOW}Warning: atom types increased from {original_n_types} to {requested_n_types}. "
            f"Remember to update Atom Types, Masses, and Pair Coeffs (and any other forcefield "
            f"sections) in the output file before running LAMMPS.{Style.RESET_ALL}"
        )

    # Show confirmation preview
    if args.mode == "block":
        confirmed = show_confirmation(
            filename=filename,
            n_chains=n_chains,
            chain_length=chain_length,
            current_labels=current_labels,
            current_lengths=current_lengths,
            requested_labels=labels,
            requested_lengths=lengths,
        )
    else:
        confirmed = show_random_confirmation(
            filename=filename,
            n_chains=n_chains,
            chain_length=chain_length,
            current_labels=current_labels,
            current_lengths=current_lengths,
            labels=labels,
            probs=probs,
            seed=args.seed,
        )

    if not confirmed:
        print("Aborted. No file written.")
        sys.exit(0)

    # Repaint and write
    if args.mode == "block":
        df_repainted = repaint(df, labels, lengths)
    else:
        df_repainted = repaint_random(df, labels, probs, seed=args.seed)

    write_output(df_repainted, lines, atoms_start, velocities_start, output_path)
    print(f"Written to: {output_path}")

    if args.mode == "random":
        hist_path = save_composition_histogram(df_repainted, labels, probs, output_path)
        print(f"Histogram:  {hist_path}")
