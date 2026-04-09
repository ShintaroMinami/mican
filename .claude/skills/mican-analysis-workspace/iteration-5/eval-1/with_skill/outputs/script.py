#\!/usr/bin/env python3
"""
MICAN protein structure comparison
Comparing: d1a05a_d1dgsa3.1.pdb vs d1a05a_d1dgsa3.2.pdb
Mode: RW (general-purpose, non-sequential alignment)
"""

from pymican import mican

# Initialize MICAN
m = mican()

# Define file paths
pdb1 = '/Users/sminami/github/mican/tests/dataset/MALISAM-ns/pdb/d1a05a_d1dgsa3.1.pdb'
pdb2 = '/Users/sminami/github/mican/tests/dataset/MALISAM-ns/pdb/d1a05a_d1dgsa3.2.pdb'

# Run alignment with RW mode (general-purpose)
print("Running MICAN alignment with RW mode (-w)...\n")
aln = m.align(pdb1=pdb1, pdb2=pdb2, options='-w')

# Calculate TM mean
tm_mean = (aln.TMscore1 + aln.TMscore2) / 2

# Standard output block
print("=" * 70)
print("ALIGNMENT RESULTS")
print("=" * 70)
print(f"TMscore  (shorter-normalized): {aln.TMscore:.4f}")
print(f"TMscore1 (pdb1-normalized):    {aln.TMscore1:.4f}  coverage: {aln.coverage1:.1f}%")
print(f"TMscore2 (pdb2-normalized):    {aln.TMscore2:.4f}  coverage: {aln.coverage2:.1f}%")
print(f"TM mean  (symmetric):          {tm_mean:.4f}")
print(f"RMSD: {aln.rmsd:.2f} Å  |  Aligned: {aln.nalign}/{aln.size1}  |  Seq ID: {aln.seq_identity:.1f}%")

# Size-difference warning
ratio = max(aln.size1, aln.size2) / min(aln.size1, aln.size2)
if ratio >= 1.5:
    print(f"\n⚠ Size ratio is {ratio:.1f}x — interpret TMscore with caution. "
          f"Check TMscore1={aln.TMscore1:.3f} vs TMscore2={aln.TMscore2:.3f}.")

# Interpretation
print("\n" + "=" * 70)
print("INTERPRETATION")
print("=" * 70)
tm = aln.TMscore
if   tm > 0.9: label = "Highly similar — typical of closely related proteins"
elif tm > 0.7: label = "Similar with notable differences — typical of distantly related proteins"
elif tm > 0.5: label = "Same fold class"
else:          label = "Low structural similarity — likely different folds"
print(f"TMscore={tm:.4f} → {label}")

# Additional details
print("\n" + "=" * 70)
print("ADDITIONAL DETAILS")
print("=" * 70)
print(f"Structure 1 size: {aln.size1} residues")
print(f"Structure 2 size: {aln.size2} residues")
print(f"Size ratio: {ratio:.2f}x")
print(f"Coverage 1: {aln.coverage1:.1f}%")
print(f"Coverage 2: {aln.coverage2:.1f}%")
