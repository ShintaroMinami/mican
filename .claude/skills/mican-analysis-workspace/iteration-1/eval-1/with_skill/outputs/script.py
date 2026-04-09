#!/usr/bin/env python3
"""Compare two PDB structures using pymican library.

This script compares two protein structures and calculates:
- TM-score
- RMSD
- Number of aligned residues
- Sequence identity
"""

import sys
from pathlib import Path

# Add mican to Python path
sys.path.insert(0, '/Users/sminami/github/mican')

from pymican import mican


def compare_structures(pdb1: str, pdb2: str) -> None:
    """Compare two PDB structures and print alignment metrics.

    Args:
        pdb1: Path to first PDB file
        pdb2: Path to second PDB file
    """
    print(f"Comparing structures:")
    print(f"  Structure 1: {pdb1}")
    print(f"  Structure 2: {pdb2}")
    print()

    # Perform alignment
    binary_path = "/Users/sminami/github/mican/build/lib.macosx-11.1-arm64-cpython-312/pymican/bin/mican"
    aligner = mican(binary=binary_path)
    alignment = aligner.align(pdb1, pdb2)

    # Get metrics
    tm_score = alignment.TMscore
    rmsd = alignment.rmsd
    n_aligned = alignment.nalign
    seq_identity = alignment.seq_identity

    # Print results
    print("Alignment Results:")
    print(f"  TM-score:           {tm_score:.4f}")
    print(f"  RMSD:               {rmsd:.4f} Å")
    print(f"  Aligned residues:   {n_aligned}")
    print(f"  Sequence identity:  {seq_identity:.2f}%")


if __name__ == "__main__":
    import os

    pdb1_path = "/Users/sminami/github/mican/tests/dataset/MALISAM-ns/pdb/d1a05a_d1dgsa3.1.pdb"
    pdb2_path = "/Users/sminami/github/mican/tests/dataset/MALISAM-ns/pdb/d1a05a_d1dgsa3.2.pdb"

    # Create output directory
    output_dir = "/Users/sminami/github/mican/.claude/skills/mican-analysis-workspace/iteration-1/eval-1/with_skill/outputs"
    os.makedirs(output_dir, exist_ok=True)

    # Redirect output to file
    output_file = os.path.join(output_dir, "result.txt")

    # Capture output
    import io
    from contextlib import redirect_stdout

    f = io.StringIO()
    with redirect_stdout(f):
        compare_structures(pdb1_path, pdb2_path)

    output = f.getvalue()

    # Write to file and print to console
    with open(output_file, 'w') as out:
        out.write(output)

    print(output)
    print(f"\nResults saved to: {output_file}")
