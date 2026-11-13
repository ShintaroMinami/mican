# MICAN Code Examples

Common setup for all snippets:

```python
import sys
import numpy as np
import pandas as pd
from pymican import mican

sys.dont_write_bytecode = True   # keep the skill directory free of __pycache__
sys.path.insert(0, "<skill base directory>/scripts")
from mican_utils import run_mican, brief_table, chain_summary, align_rank, write_superposed

m = mican()
```

## Contents
- Residue correspondence
- Complexes
- Superposition
- Batch comparison
- Sub-optimal alignments
- Preparing input files

---

## Residue correspondence

```python
CA_DIST_WARN = 4.5  # Å; aligned but not superposed beyond this


def correspondence_table(aln, residues, chain=None):
    """pdb2 counterparts of pdb1 residues (numbers as int or str)."""
    df = aln.alignment.copy()
    df['distance'] = df['distance'].astype(float)
    rows = []
    for res in residues:
        hit = df[df['residue1'] == str(res)]
        if chain is not None:
            hit = hit[hit['chain1'] == chain]
        if hit.empty:
            rows.append({'pdb1_res': str(res), 'pdb1_aa': '-', 'pdb2_res': None, 'pdb2_aa': '-',
                         'pdb2_chain': '-', 'ca_dist_A': None, 'note': 'no counterpart (unaligned)'})
            continue
        r = hit.iloc[0]
        rows.append({'pdb1_res': r['residue1'], 'pdb1_aa': r['aatype1'],
                     'pdb2_res': r['residue2'], 'pdb2_aa': r['aatype2'], 'pdb2_chain': r['chain2'],
                     'ca_dist_A': r['distance'],
                     'note': f'⚠ uncertain (>{CA_DIST_WARN} Å)' if r['distance'] > CA_DIST_WARN else 'OK'})
    return pd.DataFrame(rows)


aln = m.align('enzyme1.pdb', 'enzyme2.pdb', options='-w')
table = correspondence_table(aln, [57, 102, 195], chain='A')
print(table.to_string(index=False))
n_warn = table['note'].str.startswith('⚠').sum()
if n_warn:
    print(f"⚠ {n_warn} pair(s) are more than {CA_DIST_WARN} Å apart after superposition; "
          "treat them as uncertain correspondences.")
# table.to_csv('correspondence.csv', index=False)
```

Full table: use `aln.alignment` directly (convert `distance` to float, sort by
`residue1.astype(int)`, add the same `note` column).

---

## Complexes

```python
pdb1, pdb2 = 'complex1.cif', 'complex2.pdb'

# Whole complexes (all chains)
aln = m.align(pdb1, pdb2, options='-w')
cs = chain_summary(run_mican(pdb1, pdb2, '-w'))
print(f"C3score = {cs['C3score']:.3f}")
print(cs['chains'].to_string(index=False))   # chain1, size1, coverage1, partners {chain2: %}

# Aligned residue counts per chain pair
df = aln.alignment
print(pd.crosstab(df['chain1'], df['chain2']))

# Specific subunits: chain A of pdb1 vs chains B and C of pdb2
aln_ab = m.align(pdb1, pdb2, options='-w -c1 A -c2 BC')
```

---

## Superposition

```python
pdb1, pdb2 = 'model.cif', 'reference.pdb'

# 1. pdb1 moved into pdb2's frame, as a PDB file (all atoms, including ligands)
write_superposed(pdb1, pdb2, 'model_on_reference.pdb', options='-w')
# Open model_on_reference.pdb together with reference.pdb in PyMOL / ChimeraX.

# 2. Arbitrary coordinates of pdb1 (e.g. one ligand), moved with the same transformation
aln = m.align(pdb1, pdb2, options='-w')
xyz = np.array([[12.1, 3.4, -5.6], [13.0, 4.2, -6.1]])      # (N, 3) in pdb1's frame
xyz_on_pdb2 = aln.translate_xyz(xyz)                          # rot @ x + vec

# Ligand lines from the superposed file (works for PDB or mmCIF input)
lig = [l for l in open('model_on_reference.pdb')
       if l.startswith('HETATM') and l[17:20].strip() == 'ATP']
```

For another rank, add `-i N` to the options of `write_superposed` (it uses MICAN's own output,
which follows `-i`).

---

## Batch comparison

```python
from concurrent.futures import ThreadPoolExecutor
from itertools import combinations
from pathlib import Path

files = sorted(Path('structures/').glob('*.pdb')) + sorted(Path('structures/').glob('*.cif'))
options = '-w -f' if len(files) > 100 else '-w'


def compare(pair):
    p1, p2 = map(str, pair)
    try:
        a = m.align(p1, p2, options=options)
    except Exception:          # MICAN failed; get its message for the report
        try:
            run_mican(p1, p2, options)
            msg = 'unknown error'
        except RuntimeError as e:
            msg = str(e)
        return {'pdb1': Path(p1).name, 'pdb2': Path(p2).name, 'error': msg}
    return {'pdb1': Path(p1).name, 'pdb2': Path(p2).name, 'size1': a.size1, 'size2': a.size2,
            'TMscore1': a.TMscore1, 'TMscore2': a.TMscore2,
            'TM_mean': (a.TMscore1 + a.TMscore2) / 2, 'rmsd': a.rmsd, 'nalign': a.nalign,
            'size_ratio': max(a.size1, a.size2) / min(a.size1, a.size2), 'error': ''}


with ThreadPoolExecutor(max_workers=8) as pool:   # each call is a separate MICAN process
    df = pd.DataFrame(pool.map(compare, combinations(files, 2)))

ok = df[df['error'] == ''].sort_values('TM_mean', ascending=False)
ok['note'] = np.where(ok['size_ratio'] >= 1.5, '⚠ sizes differ: compare TMscore1 vs TMscore2', '')
print(ok.to_string(index=False))
if (df['error'] != '').any():
    print(df.loc[df['error'] != '', ['pdb1', 'pdb2', 'error']].to_string(index=False))

# Symmetric TM_mean matrix (e.g. for clustering)
names = [f.name for f in files]
mat = pd.DataFrame(np.eye(len(names)), index=names, columns=names)
for r in ok.itertuples():
    mat.loc[r.pdb1, r.pdb2] = mat.loc[r.pdb2, r.pdb1] = r.TM_mean
```

---

## Sub-optimal alignments

```python
pdb1, pdb2 = 'protein1.pdb', 'protein2.pdb'
options = '-w -n 10 -g 500'     # more candidates (-g) -> more distinct solutions

table = brief_table(run_mican(pdb1, pdb2, options))
print(table.to_string(index=False))   # rank, sTMscore, TMscore (pdb1-normalized), Dali Z, ...

# Full result for one rank (scores, alignment and rotation all from that rank)
a2 = align_rank(pdb1, pdb2, rank=2, options=options)
print(a2.TMscore1, a2.TMscore2, a2.rmsd, len(a2.alignment))

# How different is rank 2 from rank 1? Rotation angle between the two superpositions
a1 = align_rank(pdb1, pdb2, rank=1, options=options)
r = a2.translation_rot @ a1.translation_rot.T
angle = np.degrees(np.arccos(np.clip((np.trace(r) - 1) / 2, -1, 1)))
print(f"rank 2 is rotated {angle:.1f}° relative to rank 1")
```

If MICAN prints "detected only K alignments", fewer distinct solutions exist at this `-g`;
raise `-g` (e.g. 1000) to search further.

---

## Preparing input files

Keep one alternate location per atom (MICAN counts every altloc copy of a Cα as a residue):

```python
def keep_first_altloc(src, dst):
    with open(src) as fin, open(dst, 'w') as fout:
        for line in fin:
            if line.startswith(('ATOM', 'HETATM')) and line[16] not in ' A':
                continue
            fout.write(line[:16] + ' ' + line[17:] if line.startswith(('ATOM', 'HETATM')) else line)
```

For mmCIF input with an older pymican, or chain IDs longer than one character, convert or
relabel with a structure library if one is available (e.g. `gemmi convert in.cif out.pdb`).
