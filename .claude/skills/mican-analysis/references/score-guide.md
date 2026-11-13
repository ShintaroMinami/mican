# MICAN Score Guide

## Report template

```python
aln = m.align(pdb1, pdb2, options='-w')
tm_mean = (aln.TMscore1 + aln.TMscore2) / 2

print("Mode: rewiring (RW, -w)")
print(f"TM-score normalized by pdb1 ({aln.size1} res): {aln.TMscore1:.3f}   coverage {aln.coverage1:.1f}%")
print(f"TM-score normalized by pdb2 ({aln.size2} res): {aln.TMscore2:.3f}   coverage {aln.coverage2:.1f}%")
print(f"TM-score mean (symmetric summary): {tm_mean:.3f}")
print(f"RMSD {aln.rmsd:.2f} Å over {aln.nalign} aligned residues | "
      f"seq identity {aln.seq_identity:.1f}% | Dali Z {aln.DALIscore:.2f}")

ratio = max(aln.size1, aln.size2) / min(aln.size1, aln.size2)
if ratio >= 1.5:
    print(f"⚠ Sizes differ {ratio:.1f}x: the two TM-scores answer different questions (see below).")
```

## TM-score normalization

A TM-score is the alignment quality divided by the length of one of the two proteins, so the
same alignment gives two values. In pymican, `TMscore` is identical to `TMscore1` — it changes
when pdb1 and pdb2 are swapped — so always say which length a reported value uses.

| Question | Use | Why |
|----------|-----|-----|
| Are the two structures similar overall? | The lower of `TMscore1`/`TMscore2` (normalized by the larger protein) | Unaligned parts of the larger protein count against the score |
| Is the smaller structure (a domain, a fragment) contained in the larger one? | Normalized by the smaller protein | Only the smaller protein needs to be covered |
| One symmetric number (batch tables, clustering) | Mean of the two | Independent of argument order |

When sizes are similar (ratio < 1.5) the two values are close and either is fine. When they
differ, also report coverage: high TM-score with low coverage of the larger protein means a
domain-level match.

## TM-score interpretation

| TM-score | Interpretation |
|----------|----------------|
| ≥ 0.8 | Nearly identical fold; often the same protein in another state, or close homologs |
| 0.5–0.8 | Same fold |
| 0.3–0.5 | Partial similarity — a shared domain or sub-structure; check coverage |
| < 0.3 | No meaningful structural similarity (unrelated pairs typically score 0.2–0.3) |

Use these as guides, not hard cut-offs, and base the label on the normalization that fits the
question.

## Other scores

| Score | Range | Meaning |
|-------|-------|---------|
| `rmsd` | Å | Cα RMSD over aligned residues only. A small RMSD over few residues is not a good match — read it together with `nalign` and coverage |
| `nalign` | count | Number of aligned residue pairs |
| `coverage1/2` | 0–100 % | Share of each protein's residues that is aligned |
| `DALIscore` | Z-score | Dali Z; > 2 indicates significant similarity, higher is more similar |
| `SPscore` | ≥ 0 | Alignment quality score that penalizes loosely aligned pairs (higher is better) |
| `sTMscore` | 0–1 | SSE-weighted TM-score; MICAN's default criterion for ranking solutions |
| `seq_identity` | 0–100 % | Identical residues among aligned pairs (see below) |
| C3score | 0–1 | Complexes only: share of aligned residues consistent with a one-to-one chain mapping (from `chain_summary()`) |

## seq_identity and evolutionary relationship

| seq_identity | Interpretation |
|--------------|----------------|
| ≥ 30 % | Consistent with homology |
| 20–30 % | "Twilight zone" — homology is ambiguous |
| < 20 % | No detectable relationship from sequence alone |

Structure comparison cannot prove or rule out an evolutionary relationship:
- Do not conclude "analogous" or "unrelated" from TM-score or seq_identity alone; say the
  relationship is unclear.
- A high TM-score (≥ 0.8) can *suggest* homology — present it as a suggestion.
- Even ≥ 30 % identity is only consistent with homology; a definite answer needs sequence or
  phylogenetic analysis.
