# MICAN Code Examples — Basic Operations

## Residue Correspondence

### Setup

```python
from pymican import mican
import pandas as pd

CA_DIST_WARN = 4.5  # Å — pairs beyond this are structurally loose matches

m = mican()
aln = m.align('protein1.pdb', 'protein2.pdb', options='-w')

# aln.alignment columns:
#   residue1, residue2  — residue numbers (str)
#   chain1, chain2      — chain IDs
#   aatype1, aatype2    — 1-letter amino acid
#   closeness, distance — Cα distance in Å
```

### Single Residue Lookup

```python
def find_corresponding_residue(aln, residue_num, chain: str = None):
    """Return the pdb2 residue aligned to a given pdb1 residue, or None if unaligned."""
    df = aln.alignment
    mask = df['residue1'] == str(residue_num)
    if chain:
        mask &= df['chain1'] == chain
    hit = df[mask]
    if hit.empty:
        return None
    row = hit.iloc[0]
    return {'residue2': row['residue2'], 'chain2': row['chain2'],
            'aatype1': row['aatype1'], 'aatype2': row['aatype2'],
            'distance': float(row['distance'])}

result = find_corresponding_residue(aln, residue_num=42, chain='A')
if result:
    dist = result['distance']
    warning = ' ⚠ Large distance — correspondence uncertain' if dist > CA_DIST_WARN else ''
    print(f"Residue 42 (A:{result['aatype1']}) → pdb2 residue {result['residue2']} "
          f"({result['chain2']}:{result['aatype2']}), Cα dist = {dist:.2f} Å{warning}")
else:
    print("Residue 42 is outside the aligned region.")
```

### Active-Site / Key Residue Correspondence Table

```python
active_site_residues = [42, 87, 120]  # pdb1 residue numbers of interest

df = aln.alignment.copy()
df['residue1'] = df['residue1'].astype(str)
df['distance'] = df['distance'].astype(float)

rows = []
for res in active_site_residues:
    hit = df[df['residue1'] == str(res)]
    if not hit.empty:
        r = hit.iloc[0]
        rows.append({'pdb1_res': res, 'pdb1_aa': r['aatype1'], 'pdb1_chain': r['chain1'],
                     'pdb2_res': r['residue2'], 'pdb2_aa': r['aatype2'], 'pdb2_chain': r['chain2'],
                     'ca_dist_A': r['distance'],
                     'note': '⚠ uncertain (>4.5 Å)' if r['distance'] > CA_DIST_WARN else 'OK'})
    else:
        rows.append({'pdb1_res': res, 'pdb1_aa': '?', 'pdb1_chain': '-',
                     'pdb2_res': 'N/A', 'pdb2_aa': '-', 'pdb2_chain': '-',
                     'ca_dist_A': None, 'note': 'outside aligned region'})

table = pd.DataFrame(rows)
print(table.to_string(index=False))
# table.to_csv('active_site_correspondence.csv', index=False)

warned = table[table['note'].str.startswith('⚠')]
if not warned.empty:
    print(f"\n⚠ Warning: {len(warned)} residue(s) have Cα distance > {CA_DIST_WARN} Å. "
          "These pairs are structurally distant and may not represent true functional correspondence.")
```

### Full Correspondence Table

```python
df = aln.alignment.copy()
df['distance'] = df['distance'].astype(float)
df['note'] = df['distance'].apply(lambda d: '⚠ uncertain (>4.5 Å)' if d > CA_DIST_WARN else 'OK')
df = df.sort_values('residue1', key=lambda s: s.astype(int))
print(df[['chain1','residue1','aatype1','chain2','residue2','aatype2','distance','note']].to_string(index=False))
df.to_csv('full_correspondence.csv', index=False)

n_warn = (df['distance'] > CA_DIST_WARN).sum()
if n_warn:
    print(f"\n⚠ Warning: {n_warn} residue(s) have Cα distance > {CA_DIST_WARN} Å. "
          "These pairs are structurally distant and may not represent true functional correspondence.")
```

## Batch Comparison

```python
from pathlib import Path

m = mican()
pdb_files = sorted(Path('pdbs/').glob('*.pdb'))

if len(pdb_files) > 100:
    print(f"⚠ {len(pdb_files)} files found — using fast mode (-f) for screening.")

results = []
for i, p1 in enumerate(pdb_files):
    for p2 in pdb_files[i+1:]:
        aln = m.align(str(p1), str(p2), options='-w -f')
        tm_mean = (aln.TMscore1 + aln.TMscore2) / 2
        ratio = max(aln.size1, aln.size2) / min(aln.size1, aln.size2)
        results.append({
            'pdb1': p1.name, 'pdb2': p2.name,
            'TMscore': aln.TMscore, 'TMscore1': aln.TMscore1, 'TMscore2': aln.TMscore2,
            'TM_mean': tm_mean, 'rmsd': aln.rmsd, 'nalign': aln.nalign,
            'size1': aln.size1, 'size2': aln.size2,
            'size_warning': f'⚠ {ratio:.1f}x size diff' if ratio >= 1.5 else ''
        })

df = pd.DataFrame(results).sort_values('TMscore', ascending=False)
print(df.to_string(index=False))
```

## Save Alignment and Superposition Files

```python
import subprocess
from pymican import BINFILEPATH

subprocess.run([
    BINFILEPATH,
    'protein1.pdb', 'protein2.pdb',
    '-w',
    '-a', 'alignment.aln',
    '-o', 'superposed.pdb',
])
```

## Multi-Chain Structures

```python
# Align chain A of pdb1 with chain B of pdb2
aln = m.align('multi1.pdb', 'multi2.pdb', options='-w -c1 A -c2 B')
```

## Error Handling

```python
from pathlib import Path

def safe_align(m, pdb1: str, pdb2: str, options: str = '-w'):
    """Align two PDB files with basic validation."""
    for path in (pdb1, pdb2):
        if not Path(path).exists():
            raise FileNotFoundError(f"PDB file not found: {path}")
    aln = m.align(pdb1=pdb1, pdb2=pdb2, options=options)
    if aln.nalign == 0:
        print(f"⚠ No aligned residues found between {pdb1} and {pdb2}.")
        return None
    return aln
```

## Sub-optimal Alignments

> **⚠ Use `-g 500` or higher** when enumerating sub-optimal alignments — the default (`-g 50`)
> may miss valid alternative solutions. For symmetry analysis, use `-g 1000`.
>
> **⚠ For symmetry analysis**, see `symmetry-analysis.md` for how to distinguish true symmetry
> operations from false positives based on rotation angle and unit correspondence.

```python
import subprocess
from pymican import BINFILEPATH

pdb1, pdb2 = 'protein1.pdb', 'protein2.pdb'

# Top-N summary table
result = subprocess.run(
    [BINFILEPATH, pdb1, pdb2, '-w', '-n', '10', '-g', '500'],
    capture_output=True, text=True
)
in_table = False
for line in result.stdout.split('\n'):
    if 'Brief description' in line: in_table = True
    if in_table: print(line)
    if in_table and '(TMscore was' in line: break

# Residue-level detail for a specific rank (use -i RANK with -z)
def get_alignment_by_rank(pdb1, pdb2, rank, options='-w -g 500 -n 10'):
    result = subprocess.run(
        [BINFILEPATH, '-z', pdb1, pdb2] + options.split() + ['-i', str(rank)],
        capture_output=True, text=True
    )
    return [(int(l.split()[1]), int(l.split()[4]))
            for l in result.stdout.split('\n')
            if l.startswith('ALIGN') and l.split()[4] != '.']

pairs = get_alignment_by_rank(pdb1, pdb2, rank=2)
print(f"Rank 2: {len(pairs)} aligned pairs")
```
