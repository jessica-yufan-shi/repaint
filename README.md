# repaint

A CLI tool for repainting copolymer atom types in LAMMPS datafiles.

Given an existing LAMMPS molecular datafile, `repaint` reassigns atom types in each chain to follow a new block copolymer pattern — preserving all coordinates, periodic image flags, velocities, and header metadata.

## Use case

Useful when you have a simulation structure (e.g. an equilibrated diblock) and want to reuse its geometry with a different monomer sequence, such as converting a diblock into a triblock or tetrablock, without re-running a full equilibration.

## Installation

```bash
pip install -e .
```

Requires Python >= 3.10.

## Usage

```bash
repaint --file <datafile> --pattern <pattern> --lengths <lengths>
```

**Required arguments:**

| Argument | Description | Example |
|---|---|---|
| `--file` | Path to the input LAMMPS datafile | `sorted_tetrablock_data` |
| `--pattern` | Block types separated by `-` | `A-B-A` |
| `--lengths` | Block lengths separated by `-`, must sum to chain length | `10-60-10` |

**Optional arguments:**

| Argument | Description |
|---|---|
| `--output` | Output file path (default: `results/repainted_<file>_<pattern>_<lengths>`) |

**Example:**

```bash
repaint --file sorted_tetrablock_data --pattern A-B-A --lengths 10-60-10
```

Before writing any output, the tool displays a side-by-side preview of the current and requested patterns and prompts for confirmation.

## How it works

1. Reads the LAMMPS datafile and identifies chains by `mol_id`
2. Detects the existing block pattern from the first chain using run-length encoding
3. Shows a colored preview comparing the current vs. requested pattern
4. On confirmation, assigns new atom types based on the requested pattern — mapping letters to integers alphabetically (A→1, B→2, C→3, …)
5. Writes a new datafile with the updated `Atoms` section; all other sections are copied verbatim

All chains are repainted with the same pattern. Coordinates and periodic image flags are unchanged.

## Validation

`repaint` checks for the following before making any changes:

- Input file exists
- Number of blocks in `--pattern` and `--lengths` match
- Sum of `--lengths` equals the actual chain length
- If the new pattern introduces atom types not present in the original file, a warning is displayed (the LAMMPS forcefield will need to be updated separately)

## Input file format

The tool expects a standard LAMMPS molecular datafile with an `Atoms  # molecular` section:

```
atom_id  mol_id  atom_type  x  y  z  ix  iy  iz
```

A `Velocities` section is supported and preserved in the output.

## Project structure

```
repaint/
├── repaint/
│   ├── cli.py        # Argument parsing and entry point
│   ├── reader.py     # File reading and pattern detection
│   ├── painter.py    # Repainting logic and file writing
│   └── display.py    # Colored confirmation UI
├── tests/
├── examples/
│   ├── repaint.ipynb           # Notebook with manual workflow
│   ├── sorted_tetrablock_data  # Example 80-monomer tetrablock datafile
│   └── tetrablock_data         # Unsorted variant
└── pyproject.toml
```

## Dependencies

- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [colorama](https://pypi.org/project/colorama/)
