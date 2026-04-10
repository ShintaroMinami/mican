# MICAN Score Reference Guide

## Quick Start: Report All Scores

Use this pattern for all basic alignment results:

```python
aln = m.align('a.pdb', 'b.pdb', options='-w')
tm_mean = (aln.TMscore1 + aln.TMscore2) / 2

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
tm = aln.TMscore
if   tm > 0.9: label = "Highly similar structural match"
elif tm > 0.7: label = "Similar with notable structural differences"
elif tm > 0.5: label = "Same fold class"
else:          label = "Low structural similarity — likely different folds"
print(f"Interpretation: TMscore={tm:.4f} → {label}")
```

## Recommended Score by Task

| Goal | Recommended score | Reason |
|------|------------------|--------|
| Comparing proteins of similar size | `TMscore` | Reliable in this case |
| One structure is much larger (e.g., domain search) | `TMscore1` or `TMscore2` | Normalize by the structure of interest |
| Size-symmetric metric | `(TMscore1 + TMscore2) / 2` | Treats both lengths equally |
| Motif-level similarity | `sTMscore` | SSE-weighted; useful for structural motif comparison |

## TM-score Variants — Choose Carefully

MICAN outputs multiple TM-score variants normalized by different reference lengths.
**The choice matters, especially when the two structures differ significantly in size.**

| Attribute | Normalized by | Characteristic |
|-----------|--------------|----------------|
| `TMscore1` | Length of pdb1 | pdb1 perspective — tends to be higher when pdb1 is shorter |
| `TMscore2` | Length of pdb2 | pdb2 perspective — tends to be higher when pdb2 is shorter |
| `TMscore`  | Shorter structure (min) | MICAN default — can be inflated when structures differ greatly in size |

> **⚠ Size difference**: When structures differ greatly in length (e.g., 100 vs 300 residues),
> `TMscore` (normalized by the shorter) may appear high even though most of the longer structure
> is unaligned. Always check `coverage1`/`coverage2` in such cases.

## TM-score Interpretation

| TM-score | Interpretation |
|----------|----------------|
| > 0.9 | Highly similar structural match |
| > 0.7 | Similar with notable structural differences |
| > 0.5 | Same fold class. Proteins sharing only a common substructure can score in this range |
| ≤ 0.5 | Low structural similarity — likely different folds |

## Other Scores

| Score | Range | Meaning |
|-------|-------|---------|
| `rmsd` | Å | RMS Cα distance over aligned residues. Can be underestimated when `nalign` is small |
| `DALIscore` | Z-score | >2: significant similarity; >8: very similar |
| `SPscore` | 0–1 | Structural similarity accounting for alignment gaps |
| `sTMscore` | 0–1 | SSE-weighted TM-score |
| `seq_identity` | % | Sequence identity of aligned residue pairs (see interpretation below) |
| `coverage1/2` | 0–100 (%) | Fraction of each protein covered by the alignment. ⚠ Check when size ratio ≥ 1.5× |

### seq_identity and Evolutionary Relationship

| Range | Interpretation |
|-------|----------------|
| ≥ 30% | Consistent with homology (evolutionary relationship plausible) |
| < 30% | "Twilight zone" — homology is ambiguous |
| < 20% | Remote or no detectable relationship |

⚠ **Structural comparison cannot prove or disprove evolutionary relationship.**
- Never conclude "structural analog" or "no evolutionary relationship" from TM-score or seq_identity alone.
- TM-score > 0.8 may *suggest* evolutionary relationship — state as suggestion, not confirmed fact.
- seq_identity < 30%: say "relationship is unclear", not "unrelated".
- seq_identity ≥ 30%: consistent with homology, but definitive determination requires phylogenetic analysis.
