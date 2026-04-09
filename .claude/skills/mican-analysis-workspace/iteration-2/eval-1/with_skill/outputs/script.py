#!/usr/bin/env python3
"""Compare two PDB structures using pymican library with RW (Rewiring) mode.

This script follows the MICAN analysis skill instructions and compares:
- tests/dataset/MALISAM-ns/pdb/d1a05a_d1dgsa3.1.pdb
- tests/dataset/MALISAM-ns/pdb/d1a05a_d1dgsa3.2.pdb

Mode: RW (-w) - non-sequential alignment for general-purpose comparison
"""

from pathlib import Path
from pymican import mican


def main() -> None:
    """Compare two PDB structures and display detailed alignment metrics."""
    # Define paths
    pdb1 = '/Users/sminami/github/mican/tests/dataset/MALISAM-ns/pdb/d1a05a_d1dgsa3.1.pdb'
    pdb2 = '/Users/sminami/github/mican/tests/dataset/MALISAM-ns/pdb/d1a05a_d1dgsa3.2.pdb'

    # Display header
    print(f"Comparing {Path(pdb1).name} and {Path(pdb2).name}")
    print(f"Mode: RW (Rewiring) - non-sequential alignment")
    print("=" * 70)

    # Instantiate MICAN and run alignment
    m = mican()
    aln = m.align(pdb1=pdb1, pdb2=pdb2, options='-w')

    # Display key scores
    print(f"\nAlignment Results:")
    print(f"TM-score:        {aln.TMscore:.4f}  (>0.5 = similar fold)")
    print(f"TM-score1:       {aln.TMscore1:.4f} (normalized by protein 1 length)")
    print(f"TM-score2:       {aln.TMscore2:.4f} (normalized by protein 2 length)")
    print(f"RMSD:            {aln.rmsd:.2f} Å")
    print(f"Aligned:         {aln.nalign} residues")
    print(f"Protein 1 size:  {aln.size1} residues")
    print(f"Protein 2 size:  {aln.size2} residues")
    print(f"Coverage 1:      {aln.coverage1:.1%}")
    print(f"Coverage 2:      {aln.coverage2:.1%}")
    print(f"Seq Identity:    {aln.seq_identity:.1f}%")
    print(f"DALI score:      {aln.DALIscore:.2f}")
    print(f"SP score:        {aln.SPscore:.4f}")

    # Display first 10 aligned residue pairs
    print(f"\nFirst 10 aligned residue pairs:")
    print(aln.alignment.head(10).to_string(index=False))

    # Interpretation
    print(f"\n{'=' * 70}")
    print("Interpretation:")
    if aln.TMscore > 0.5:
        print("These proteins have SIMILAR FOLDS (TM-score > 0.5)")
    elif aln.TMscore > 0.2:
        print("These proteins have SOME STRUCTURAL SIMILARITY (0.2 < TM-score < 0.5)")
    else:
        print("These proteins are STRUCTURALLY UNRELATED (TM-score < 0.2)")


if __name__ == "__main__":
    main()
