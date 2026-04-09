---
name: mican-analysis
description: >
  Use this skill for any protein 3D structure analysis using MICAN: alignment, TM-score,
  RMSD, structural similarity, superposition, batch comparison, residue correspondence,
  conformational change (apo/holo), internal symmetry detection, rotational symmetry (Cn),
  domain detection, sub-optimal alignments, cyclic permutation, tandem repeat.
  Trigger on phrases like "compare PDB files", "align proteins", "calculate TM-score",
  "run MICAN", "structural similarity", "does this protein have symmetry", "analyze
  conformational change", "show alternative alignments", or whenever PDB file paths are
  provided. Symmetry analysis, conformational change quantification, and domain detection
  are all within scope — trigger even when the task goes beyond simple alignment.
---

# MICAN Protein Structure Analysis Skill

MICAN compares protein 3D structures via structure alignment. This skill contains the
workflows needed to execute each type of analysis correctly.

**Code templates** are in the reference files — read them when you need exact code:
- `references/score-guide.md` — TM-score interpretation and standard output block
- `references/code-examples.md` — residue correspondence, batch, superposition, sub-optimal
- `references/symmetry-analysis.md` — self-comparison, rotation matrix, cyclic permutation, domain detection
- `references/conformational-analysis.md` — per-residue Cα displacement

---

## Setup

`mican` is a class. Always instantiate it first:

```python
from pymican import mican
m = mican()
aln = m.align(pdb1='a.pdb', pdb2='b.pdb', options='-w')
```

`aln` has these key attributes:
- `TMscore`, `TMscore1`, `TMscore2` — TM-scores (normalized by min/pdb1/pdb2 length)
- `coverage1`, `coverage2` — alignment coverage in % (0–100, **not** 0–1)
- `rmsd`, `nalign`, `size1`, `size2`, `seq_identity`
- `alignment` — pandas DataFrame; columns: `residue1`, `residue2` (str), `chain1`, `chain2`, `aatype1`, `aatype2`, `distance` (Cα in Å)

---

## Workflow A: Compare Two Proteins

**Step 1 — Choose alignment mode** (if not specified, use RW and announce it):

| Situation | Mode | Option |
|-----------|------|--------|
| Same protein, different states | Sequential (SQ) | `-s` |
| Unknown relationship (default) | General-purpose (RW) | `-w` |
| Want reverse/antiparallel SSEs too | Exhaustive (RR) | `-r` |
| Large-scale screening | add Fast | `-f` |

**Step 2 — Run and report all TM-score variants** (use Standard Output Block in `score-guide.md`):

```python
aln = m.align(pdb1='a.pdb', pdb2='b.pdb', options='-w')
tm_mean = (aln.TMscore1 + aln.TMscore2) / 2
# Print TMscore, TMscore1+coverage1, TMscore2+coverage2, TM mean, RMSD, nalign, seq_identity
# Size warning if max(size1,size2)/min(size1,size2) >= 1.5
```

**Step 3 — Interpret results:**

| TM-score | Meaning |
|----------|---------|
| > 0.9 | Highly similar |
| > 0.7 | Similar, notable differences (distantly related) |
| > 0.5 | Same fold class |
| ≤ 0.5 | Low similarity — likely different folds |

- **Size ratio ≥ 1.5**: TMscore (shorter-normalized) is inflated. Use TMscore1 or TMscore2 instead, or TM mean for a symmetric view.
- **Low coverage**: only part of one protein aligns — possible domain match, not whole-protein similarity.

---

## Workflow B: Residue Correspondence / Active-Site Mapping

**Step 1** — Run alignment (RW recommended).

**Step 2** — Query `aln.alignment`. Key trap: `residue1` and `residue2` are **strings** — use `str(N)` for lookup.

**Step 3** — Flag unreliable pairs: Cα distance > 4.5 Å means the aligned residues are spatially distant and should not be treated as functionally corresponding. Add a `note` column: `'⚠ uncertain (>4.5 Å)'` or `'OK'`.

**Step 4** — Print table and warn if any pairs exceed 4.5 Å. Save to CSV if requested.

Use code from `references/code-examples.md` — Residue Correspondence.

---

## Workflow C: Internal Symmetry and Rotational Symmetry

This is the most complex workflow. Follow every step in order.

### Step 1: Self-comparison (always start here)

Run the protein against itself with large `-g` and `-n`:

```python
# Use subprocess, not pymican.mican() — need raw output including rotation matrices
result = subprocess.run(
    [BINFILEPATH, pdb, pdb, '-w', '-g', '1000', '-n', '20'],
    capture_output=True, text=True
)
```

Read the "Brief description" table from stdout. **DALI Z > 2 in Rank 2+** means statistically significant internal symmetry.

### Step 2: Classify sub-optimal solutions — True vs False Positives

Before any further analysis, separate the solutions:

**True symmetry operations** (analyze these):
- Rotation angle **> 30°**
- Residue correspondence between **different** repeat units (u0→u1, u0→u2, etc.)
- Solutions appear in **pairs** (forward + reverse direction)

**False positives** (discard these):
- Rotation angle **< 20°** — nearly identity transformation; MICAN found spatially nearby regions, not a symmetry operation
- Correspondence **within the same unit** (u0→u0)

To get rotation angle for each rank, use `-z -i RANK` and parse rotation matrix lines starting with ` 1   `, ` 2   `, ` 3   `. Then: `cos θ = (trace(R) − 1) / 2`.

### Step 3: Rotation matrix analysis — identify Cn symmetry

For each true-symmetry rank:
1. **Shared axis?** If all ranks have rotation axes within ~1° of each other → single true symmetry axis exists.
2. **Cn match?** Check if rotation angles equal `k/n × 360°` (within ~1°) for some integer n.
3. **Pure rotation or screw?** Compute `t∥ = t · axis`. If ≈ 0 → pure Cn rotation. If ≠ 0 → screw symmetry.
4. **Pseudo higher-order symmetry?** If some ranks show angles that are close to (but not matching within ~1°) a higher-order Cn, note this as possible pseudo Cn character. Do not assert a specific fold class (e.g., β-propeller) from rotation angles alone — that requires additional evidence beyond MICAN's structural comparison.

Use code from `references/symmetry-analysis.md` §2.

### Step 4: Cyclic permutation — algebraic proof of Cn

Map each residue to its repeat unit index. Check whether the unit-level mapping across ranks forms all cyclic shifts {0, 1, ..., n−1}. If yes → Cn is algebraically proven.

Estimate repeat unit size from GCD of alignment shifts across ranks.

Use code from `references/symmetry-analysis.md` §3.

### Step 5: Check for multiple domains first (do this before Step 3 if uncertain)

Signs that the protein has multiple structural domains:
- High-scoring ranks (Rank 2–K) cover only **part** of the chain (res1 max ≤ X)
- Lower ranks suddenly expand to the full chain (res1 max jumps)
- Some ranks show angles that don't fit any Cn hypothesis → those ranks span a domain boundary

If domain signals appear: use `extract_domain()` to split the PDB at the estimated boundary, then analyze each domain separately (restart from Step 1 for each domain).

Use code from `references/symmetry-analysis.md` §4.

---

## Workflow D: Conformational Change (Apo/Holo)

**Step 1** — Align on a tight core using `-q` threshold:
- Use **SQ mode** (`-s`) — same protein, sequence order is conserved
- `-q 3.0` is a good default (2.0 for very tight core, 4.5 for loose)
- The core = residues that don't move between states (the anchor)

**Step 2** — Apply the rotation matrix to ALL residues (not just aligned ones), then measure Cα displacement:
- `in_core=True` residues: small displacement — structurally invariant anchor
- `in_core=False` with large displacement (> 5 Å): moved regions — conformational change
- Large jump between core RMSD (~2 Å) and non-core RMSD (>> 5 Å) → rigid-body domain movement

Use `analyze_conformational_change()` from `references/conformational-analysis.md`.

---

## Key CLI Options

| Option | Default | Notes |
|--------|---------|-------|
| `-n INT` | 5 | Solutions to return; use 20 for symmetry analysis |
| `-g INT` | 50 | GH candidates; use 1000 for symmetry/sub-optimal analysis |
| `-q FLOAT` | — | Cα distance cutoff (tighter = more precise core) |
| `-i INT` | — | Retrieve rank N solution (requires `-z`) |
| `-z` | off | Machine-readable output; **required** for `-i` and rotation matrix parsing |
| `-c1 ID` / `-c2 ID` | — | Chain ID for pdb1/pdb2 |
| `-l INT` | 3 | Minimum segment length |
