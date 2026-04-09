---
name: mican-analysis
description: >
  Use this skill whenever the user wants to compare, align, or analyze protein 3D structures
  using MICAN. Trigger when the user mentions: comparing PDB files, protein structure alignment,
  TM-score, RMSD, non-sequential alignment, structural similarity, superposition of proteins,
  batch structure comparison, or asks things like "compare these two proteins", "how similar are
  these structures", "align protein A and B", "calculate TM-score", or "run MICAN on these files".
  Also trigger when the user provides PDB file paths and asks anything about structural analysis.
---

# MICAN Protein Structure Analysis Skill

MICAN is a protein structure alignment tool that handles non-sequential, multi-chain,
and reverse-direction alignments. Use this skill to run MICAN via Python and interpret results.

## Reference Files

Read these when you need details beyond what is in this file:
- `references/score-guide.md` — TM-score variants, interpretation table, recommended score by task, standard output block
- `references/code-examples.md` — ready-to-use code for correspondence tables, batch comparison, superposition, multi-chain, error handling

## Quick Patterns

**"Compare these two PDB files"** (no mode specified) → Announce RW default and proceed
**"Which mode should I use?"** → Explain using plain-language descriptions (see Alignment Modes), never bare option names
**"Align many proteins"** → See `references/code-examples.md` — Batch Comparison
**"Show me which residues aligned"** → See `references/code-examples.md` — Residue Correspondence
**"Save the superposed structure"** → See `references/code-examples.md` — Save Alignment and Superposition Files
**"The proteins have a circular permutation"** → Use general-purpose (`-w`) or exhaustive (`-r`) alignment
**"Which residue in protein B corresponds to active-site residue N?"** → See `references/code-examples.md` — Single Residue Lookup
**"Generate a correspondence table"** → See `references/code-examples.md` — Active-Site / Key Residue Correspondence Table
**"Show the full residue mapping"** → See `references/code-examples.md` — Full Correspondence Table
**"Interpret the TM-score"** → See `references/score-guide.md`

## Setup

Always import and instantiate as follows — `mican` is a class, not a standalone function:

```python
from pymican import mican   # import the mican class
m = mican()                 # instantiate (uses bundled binary automatically)
aln = m.align(pdb1='protein1.pdb', pdb2='protein2.pdb', options='-w')
```

If `pymican` is not installed, build and install it first:
```bash
make  # builds pymican/bin/mican from src/
python setup.py install
```

## Alignment Mode (Step 0)

**If the user has not specified a mode, default to RW and announce it before running:**

> "No alignment mode was specified — running with the general-purpose non-sequential alignment (RW mode)."

Then proceed immediately without asking for confirmation.

**If the user asks which mode to use**, present the options:

| Mode | Description |
|------|-------------|
| **Sequential alignment** (SQ, `-s`) | Matches residues in order. Best for homologous proteins with conserved sequence order. |
| **General-purpose alignment** (RW, `-w`) | Finds structurally similar regions regardless of sequence order. Standard choice when relationship is unknown. **Default.** |
| **Exhaustive alignment** (RR, `-r`) | Like RW, but also explores reverse/antiparallel SSE pairings. More thorough but slower. |
| Fast mode (`-f`) | Add to any mode for large-scale screening. |

Always use plain-language descriptions (not raw option names like `-w`) when explaining modes to the user. Option names may appear in code only.

## Basic Alignment

```python
from pymican import mican

m = mican()
aln = m.align(pdb1='protein1.pdb', pdb2='protein2.pdb', options='-w')
```

After alignment, always report all TM-score variants and check for size differences.
Read `references/score-guide.md` — Standard Output Block for the complete output pattern.

Key attributes:
- `aln.TMscore` / `aln.TMscore1` / `aln.TMscore2` — TM-score (see score-guide.md for which to use)
- `aln.coverage1` / `aln.coverage2` — alignment coverage in % (0–100)
- `aln.rmsd` — RMS Cα distance over aligned residues (Å)
- `aln.nalign` / `aln.size1` / `aln.size2` — aligned and total residue counts
- `aln.seq_identity` — sequence identity of aligned pairs (%)
- `aln.alignment` — pandas DataFrame of aligned residue pairs (see Residue Correspondence below)

## TM-score — Quick Decision

| Situation | Use |
|-----------|-----|
| Structures similar in size (ratio < 1.5) | `TMscore` |
| One structure much larger | `TMscore1` or `TMscore2` (normalize by structure of interest) |
| Size-symmetric comparison | `(TMscore1 + TMscore2) / 2` |

**Always check `coverage1`/`coverage2` when sizes differ significantly.**
For full interpretation table, see `references/score-guide.md`.

## Residue Correspondence

`aln.alignment` is a pandas DataFrame — each row is one aligned residue pair.

| Column | Description |
|--------|-------------|
| `residue1`, `residue2` | Residue numbers (str) |
| `chain1`, `chain2` | Chain IDs |
| `aatype1`, `aatype2` | 1-letter amino acid |
| `distance` | Cα distance in Å |

**Distance reliability threshold: `CA_DIST_WARN = 4.5` Å**
Pairs with Cα distance > 4.5 Å are structurally loose — always flag them with `⚠` in output
and print a summary warning so the user is not misled about functional correspondence.

For complete code patterns (single lookup, key-residue table, full table), read
`references/code-examples.md` — Residue Correspondence.

## Key CLI Options

| Option | Description |
|--------|-------------|
| `-n INT` | Number of solutions to return (default: 5) |
| `-g INT` | GH candidates (default: 50; increase for better coverage) |
| `-l INT` | Minimum segment length (default: 3) |
| `-c1 ID` / `-c2 ID` | Chain ID for pdb1/pdb2 |
| `-z` | Machine-readable output (used internally by pymican) |
