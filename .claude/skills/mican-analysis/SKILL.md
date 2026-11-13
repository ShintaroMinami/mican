---
name: mican-analysis
description: >
  Compare protein 3D structures with MICAN (pymican): structure alignment, TM-score, RMSD,
  Dali Z, residue-to-residue correspondence (e.g. mapping active-site residues between
  proteins), multi-chain complexes (chain selection, chain-to-chain correspondence),
  superposed coordinates, batch / all-vs-all comparison, and alternative (sub-optimal)
  alignments. Handles non-sequential (rewired, circularly permuted) and reverse-direction
  alignments, and PDB or mmCIF input. Use this skill whenever the user wants to compare,
  align or superpose protein structures, asks how similar two structures are, wants the
  residues of one structure that correspond to residues of another, or gives .pdb/.cif files
  to compare — and whenever MICAN or pymican is mentioned, even if the tool is not named.
---

# MICAN Structure Comparison

MICAN aligns two protein structures and scores their similarity with TM-score and related
measures. Unlike sequential aligners, it can also match secondary-structure segments that
appear in a different order (circular permutations, swapped elements) or run in the opposite
direction, and it aligns whole multi-chain complexes at once.

Work through **pymican** (the Python wrapper). For what pymican does not expose — MICAN's error
messages, ranks other than the best, per-chain results, and superposed coordinate files — use
the helpers in `scripts/mican_utils.py` (in this skill's base directory) instead of parsing
MICAN's output yourself:

```python
import sys
sys.dont_write_bytecode = True   # keep the skill directory free of __pycache__
sys.path.insert(0, "<skill base directory>/scripts")
from mican_utils import run_mican, brief_table, chain_summary, align_rank, write_superposed
```

Read these when needed:
- `references/score-guide.md` — before reporting or interpreting scores (report template,
  TM-score normalization, what seq_identity can and cannot tell you).
- `references/code-examples.md` — code for residue correspondence, complexes, superposition,
  batch comparison and sub-optimal alignments.

## Setup

If `import pymican` fails, run `python <skill base directory>/scripts/setup_env.py`. It installs
pymican with pip; building the MICAN binary needs a C compiler (gcc or clang).

## pymican basics

```python
from pymican import mican
m = mican()
aln = m.align('a.pdb', 'b.cif', options='-w')   # file paths, or file contents as strings
```

| Attribute | Meaning |
|-----------|---------|
| `TMscore1`, `TMscore2` | TM-score normalized by the length of pdb1 / pdb2 |
| `TMscore` | Same value as `TMscore1` (normalized by **pdb1**, not by the shorter chain) |
| `coverage1`, `coverage2` | % of pdb1 / pdb2 residues aligned (0–100, not 0–1) |
| `rmsd`, `nalign` | Cα RMSD (Å) over the `nalign` aligned residues |
| `size1`, `size2` | Number of residues (Cα atoms) read from each file |
| `seq_identity` | % identical residues among aligned pairs |
| `DALIscore`, `SPscore`, `sTMscore` | Dali Z-score, SP-score, SSE-weighted TM-score |
| `alignment` | DataFrame: `residue1`, `residue2`, `chain1`, `chain2`, `aatype1`, `aatype2`, `distance` |
| `translation_rot`, `translation_vec`, `translate_xyz()` | Transformation that moves pdb1 onto pdb2 |

Things that commonly go wrong:
- **All `alignment` columns are strings.** Compare residues with `str(42)` and convert
  `distance` with `.astype(float)` before doing arithmetic.
- **pymican hides MICAN's errors.** A failed run surfaces as a `KeyError` (e.g. `'mode'`) or a
  parsing error. Call `run_mican(pdb1, pdb2, options)` to see MICAN's actual message (missing
  file, no Cα atoms, unsupported chain ID, ...).
- **`-i N` mixes ranks in pymican.** With `options='-w -i 2'`, the alignment, `TMscore1/2`,
  coverage and rotation are for rank 2, but `TMscore`, `rmsd`, `nalign`, `seq_identity`,
  `DALIscore`, `SPscore` and `sTMscore` stay at rank 1. Use `align_rank()` for any rank other
  than 1.

## Choose the alignment mode

If the user does not specify one, use RW (`-w`) and say so: it finds sequential alignments too,
and also catches rearranged topologies that a sequential search would miss.

| Situation | Mode | Option |
|-----------|------|--------|
| General comparison, unknown relationship (default) | Rewiring (RW) | `-w` |
| Sequence order must be preserved (same protein in different states, comparing with TM-align-style results) | Sequential (SQ) | `-s` |
| Also allow segments in the reverse direction (most exhaustive, slower) | Rewiring + reverse (RR) | `-r` |
| Large screening where speed matters | add fast mode | `-f` |

## Workflow A: Compare two structures

1. Run the alignment in the chosen mode.
2. Report both normalizations — `TMscore1` with `coverage1` and `TMscore2` with `coverage2` —
   plus RMSD, `nalign`, `seq_identity` and Dali Z, using the template in `score-guide.md`.
   A single TM-score is ambiguous because it depends on which length it was divided by.
3. Interpret (details in `score-guide.md`):

   | TM-score | Meaning |
   |----------|---------|
   | ≥ 0.8 | Nearly identical fold; often the same protein or close homologs |
   | 0.5–0.8 | Same fold |
   | 0.3–0.5 | Partial similarity (a shared domain or sub-structure) |
   | < 0.3 | No meaningful similarity (unrelated pairs typically score 0.2–0.3) |

   - When sizes differ (ratio ≥ 1.5), say which normalization answers the question: normalizing
     by the smaller protein asks "is it contained in the larger one?"; normalizing by the larger
     asks "are they similar overall?". Low coverage of the larger protein means a domain-level
     match, not whole-protein similarity.
   - Structure similarity alone does not establish or exclude an evolutionary relationship; see
     `score-guide.md` before commenting on homology.

## Workflow B: Residue correspondence / active-site mapping

1. Align (RW unless sequence order matters).
2. Look up residues in `aln.alignment` (strings!). Residues missing from the table are
   unaligned — report them as "no counterpart", not as an error.
3. Flag pairs whose Cα `distance` exceeds 4.5 Å as uncertain: they were placed in the
   alignment but are not spatially superposed, so they should not be treated as functionally
   equivalent.
4. Print a table (pdb1 residue, aa, chain → pdb2 residue, aa, chain, distance, note) and save
   CSV if asked. Code: `code-examples.md` § Residue correspondence.

## Workflow C: Multi-chain complexes

By default MICAN reads **all chains** of each file and aligns the complexes as wholes; one chain
of pdb1 may be matched to parts of several chains of pdb2.

1. Align the complexes; report the usual scores plus **C3score** and the per-chain results from
   `chain_summary(run_mican(...))`.
   - C3score is the fraction of aligned residues that follow a one-to-one chain mapping. 1.0
     means every chain of pdb1 corresponds to a single chain of pdb2; low values mean the
     alignment jumps between chains, so the subunit correspondence is unclear.
   - Per-chain rows give, for each pdb1 chain, its coverage and how it splits over pdb2 chains.
2. To compare specific subunits, select chains with `-c1` / `-c2` (`-c1 A`, `-c1 AB` for
   chains A and B). The `chain1` / `chain2` columns of `aln.alignment` tell which chain each
   aligned residue belongs to.

Code: `code-examples.md` § Complexes.

## Workflow D: Superposed coordinates

The transformation moves **pdb1 onto pdb2**; pdb2's coordinates never change.
- Structure file for a viewer: `write_superposed(pdb1, pdb2, 'pdb1_on_pdb2.pdb', '-w')` writes
  pdb1's atoms in pdb2's frame; open it together with the original pdb2.
- Arbitrary coordinates (e.g. a ligand or a site from pdb1): `aln.translate_xyz(xyz)` with an
  (N, 3) array.
- MICAN's own `-o file` is a RasMol script (MODEL 1 = moved pdb1, MODEL 2 = pdb2) and `-a file`
  writes the alignment as text; use these only when the user asks for MICAN's file formats.

Code: `code-examples.md` § Superposition.

## Workflow E: Batch comparison

Loop over pairs with pymican; each call is a separate process, so a thread pool parallelizes
well. Keep `TMscore1` and `TMscore2` for every pair (the score is asymmetric), flag pairs whose
size ratio is ≥ 1.5, catch failures per pair instead of aborting the batch, and add `-f` for
more than ~100 structures. Sort or cluster by the normalization that fits the question (the
mean of both is a good symmetric default). Code: `code-examples.md` § Batch comparison.

## Workflow F: Alternative (sub-optimal) alignments

MICAN ranks several distinct solutions (different superpositions), not just the best one. They
are useful when the top alignment is ambiguous — repeats, multi-domain proteins, or when the user
asks for other possible matches.

- `-n N` sets how many solutions to report (default 5); `-g K` sets how many candidates are
  explored (default 50). Raise `-g` to 200–1000 when fewer than N solutions come back.
- `brief_table(run_mican(pdb1, pdb2, '-w -n 10 -g 500'))` gives the scores of every rank;
  `align_rank(pdb1, pdb2, rank=2, options='-w -n 10 -g 500')` gives the full result of one rank.

Code: `code-examples.md` § Sub-optimal alignments.

## Input files

- PDB and PDBx/mmCIF are both accepted; the format is detected from the file contents.
  mmCIF needs pymican ≥ 1.13.0 — with older versions MICAN reports "No CA atoms" for a .cif
  file; then upgrade (`pip install -U pymican`) or convert the file to PDB.
- Only the first model is read, and each Cα atom (ATOM or HETATM, e.g. MSE) counts as a residue.
  Alternate locations of a Cα are each counted, which inflates `size1`/`size2`; keep a single
  altloc first if the residue count looks too large.
- mmCIF files use author chain IDs. Chain IDs longer than one character are not supported: select
  single-character chains with `-c1`/`-c2`, or the run stops with an error.

## Key CLI options

| Option | Default | Notes |
|--------|---------|-------|
| `-s` / `-w` / `-r` | SQ | Alignment mode (see above) |
| `-f` | off | Fast mode for screening (same as `-g 15`) |
| `-c1 IDs` / `-c2 IDs` | all chains | Chains of pdb1 / pdb2, e.g. `-c1 AB` |
| `-n INT` | 5 | Number of solutions reported |
| `-g INT` | 50 | Candidates explored; raise for more or better sub-optimal solutions |
| `-i INT` | 1 | Report the i-th solution in detail (use `align_rank()` from Python) |
| `-q REAL` | 10.0 | Maximum Cα distance (Å) for aligned pairs; lower it for a stricter core |
| `-l INT` | 3 | Minimum segment length |
| `-t INT` | 0 | Score used to rank solutions: 0 sTMscore, 1 TMscore, 2 Dali Z |
| `-o` / `-a` / `-m` file | — | Superposition (RasMol script) / alignment / matrix files |
