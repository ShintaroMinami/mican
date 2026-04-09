# MICAN Conformational Change Analysis

Quantify how much each residue moves between two states (e.g., apo vs holo) by:
1. Aligning on a structurally conserved core (`-q` threshold)
2. Applying the rotation matrix to all residues
3. Measuring per-residue Cα displacement

## Usage

```python
import subprocess
import numpy as np
import pandas as pd
from pymican import BINFILEPATH

def analyze_conformational_change(pdb_state1: str, pdb_state2: str,
                                   q_threshold: float = 3.0) -> pd.DataFrame:
    """
    Superpose state1 onto state2 using a tight core alignment (-q threshold),
    then compute per-residue Cα displacement for all common residues.

    Returns DataFrame with columns: residue, ca_dist_A, in_core
    Tip: use -s (sequential) mode — states of the same protein share sequence order.
    """
    # Step 1: core alignment and rotation matrix
    result = subprocess.run(
        [BINFILEPATH, '-z', pdb_state1, pdb_state2,
         '-s', f'-q{q_threshold}', '-n', '1', '-g', '100'],
        capture_output=True, text=True
    )
    core_residues = set()
    rot = np.zeros((3, 3)); vec = np.zeros(3)
    for line in result.stdout.split('\n'):
        if line.startswith('ALIGN') and line.split()[4] != '.':
            core_residues.add(int(line.split()[1]))
        if line.startswith(' 1   '):
            _, vec[0], rot[0,0], rot[0,1], rot[0,2] = map(float, line.split())
        if line.startswith(' 2   '):
            _, vec[1], rot[1,0], rot[1,1], rot[1,2] = map(float, line.split())
        if line.startswith(' 3   '):
            _, vec[2], rot[2,0], rot[2,1], rot[2,2] = map(float, line.split())

    print(f"Core residues (q={q_threshold}Å): {len(core_residues)}")

    # Step 2: read Cα coordinates
    def read_ca(path):
        coords = {}
        with open(path) as f:
            for line in f:
                if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                    resnum = int(line[22:26])
                    coords[resnum] = np.array([float(line[30:38]),
                                               float(line[38:46]),
                                               float(line[46:54])])
        return coords

    ca1 = read_ca(pdb_state1)
    ca2 = read_ca(pdb_state2)

    # Step 3: apply rotation, compute displacement
    rows = []
    for res in sorted(set(ca1) & set(ca2)):
        xyz_sup = rot @ ca1[res] + vec
        dist = float(np.linalg.norm(xyz_sup - ca2[res]))
        rows.append({'residue': res, 'ca_dist_A': dist, 'in_core': res in core_residues})

    df = pd.DataFrame(rows)

    # Step 4: summary
    core_df    = df[df['in_core']]
    noncore_df = df[~df['in_core']]
    print(f"Core    RMSD: {np.sqrt((core_df['ca_dist_A']**2).mean()):.2f} Å")
    print(f"Non-core RMSD: {np.sqrt((noncore_df['ca_dist_A']**2).mean()):.2f} Å  "
          f"max: {noncore_df['ca_dist_A'].max():.2f} Å")

    mobile = noncore_df[noncore_df['ca_dist_A'] > 5.0].sort_values('ca_dist_A', ascending=False)
    if not mobile.empty:
        print(f"Residues displaced > 5 Å (n={len(mobile)}): "
              f"{mobile['residue'].tolist()[:10]}{'...' if len(mobile) > 10 else ''}")
    return df

# Usage
df = analyze_conformational_change('apo.pdb', 'holo.pdb', q_threshold=3.0)
# df columns: residue, ca_dist_A, in_core
```

## Interpreting results

- **`in_core=True` residues** — structurally conserved anchor; changes in these are minor
- **`in_core=False` with large displacement** — regions that move upon conformational change
- A large jump in displacement between core and non-core (e.g., core RMSD 2 Å vs non-core RMSD 20+ Å) indicates a rigid-body domain movement

## Choosing `-q` threshold

| `-q` value | Effect |
|-----------|--------|
| 2.0 Å | Very tight — only the most invariant core |
| 3.0 Å | Recommended default — good balance |
| 4.5 Å | Loose — larger core but less reliable anchor |

When `-q` is too loose, core RMSD will be inflated and displacement contrast is reduced.

## Visualization (matplotlib)

```python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 1, figsize=(14, 7), gridspec_kw={'height_ratios': [3, 1]})
colors = ['steelblue' if c else 'tomato' for c in df['in_core']]
axes[0].bar(df['residue'], df['ca_dist_A'], color=colors, width=1.0)
axes[0].axhline(5.0,  color='orange', linestyle='--', linewidth=1.2, label='5 Å')
axes[0].axhline(10.0, color='red',    linestyle='--', linewidth=1.2, label='10 Å')
axes[0].set_xlabel('Residue number'); axes[0].set_ylabel('Cα displacement (Å)')
axes[0].legend()
for _, row in df.iterrows():
    axes[1].axvspan(row['residue']-0.5, row['residue']+0.5,
                    color='steelblue' if row['in_core'] else 'tomato', alpha=0.7)
axes[1].set_yticks([]); axes[1].set_xlabel('Residue')
plt.tight_layout()
plt.savefig('conformational_change.png', dpi=150, bbox_inches='tight')
```
