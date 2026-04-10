# MICAN Conformational Change Analysis

Quantify how much each residue moves between two states (e.g., apo vs holo) by:
1. Aligning on a structurally conserved core (`-q` threshold)
2. Applying the rotation matrix to all residues
3. Measuring per-residue Cα displacement

## Step 0: Assess Protein Type Before Choosing Strategy

**Before running the analysis, check whether the protein has internal repeats or symmetry.**
The appropriate strategy depends on the structural type:

| Protein type | Strategy |
|-------------|----------|
| **Non-repeat protein** | Analyze whole chain directly → proceed to Step 1 then Step 2 |
| **Repeat / symmetric protein** | First run symmetry analysis (see `symmetry-analysis.md`), identify repeat units, extract each unit, then analyze each unit independently → Skip to Step 3 |

Why this matters for repeat proteins:
- A tandem-repeat or symmetric protein undergoing conformational change may shift the
  register of repeats between states (e.g., one state is C4, another is slightly distorted).
  Applying whole-chain alignment will mix inter-repeat displacement with intra-repeat
  displacement, obscuring the true conformational change.
- Extract each repeat unit from both states and compare unit-by-unit.

## Step 1: Choose `-q` Threshold

For **non-repeat proteins**, use a tight threshold — it defines the "rigid core" that anchors
the superposition, and a loose threshold dilutes the signal:

| `-q` value | When to use |
|-----------|-------------|
| **~3.0 Å** | **Recommended for most conformational change analyses.** Captures the structurally conserved core while excluding clearly mobile regions. |
| 2.0 Å | Stricter — use when the conformational change is small and you need a tighter anchor |
| > 4.5 Å | Too loose for most cases — core RMSD becomes inflated and the anchor is unreliable |

> **Rule of thumb:** Start with `-q 3.0`. If too few core residues are found (< 10% of chain),
> relax to 4.0. If the core RMSD seems suspiciously large, tighten to 2.0.

## Step 2: Run Per-Residue Displacement Analysis

```python
import subprocess
import numpy as np
import pandas as pd
from pymican import BINFILEPATH

def analyze_conformational_change(pdb_state1: str, pdb_state2: str,
                                   q_threshold: float = 3.0) -> pd.DataFrame:
    """
    Superpose state1 onto state2 using a conserved core alignment (-q threshold),
    then compute per-residue Cα displacement for all common residues.

    Returns DataFrame with columns: residue, ca_dist_A, in_core

    Notes:
    - Use -s (sequential) mode — states of the same protein share sequence order.
    - Default q_threshold=3.0 Å is recommended for general conformational change analysis.
    - For repeat/symmetric proteins, extract individual repeat units first and call
      this function on each unit separately (see Step 3 in this file).
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

    n_core = len(core_residues)
    print(f"Core residues (q={q_threshold}Å): {n_core}")
    if n_core == 0:
        print(f"⚠ No core residues found at q={q_threshold} Å. Try relaxing to q=5.0.")
        return pd.DataFrame(columns=['residue', 'ca_dist_A', 'in_core'])

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
    if not noncore_df.empty:
        print(f"Non-core RMSD: {np.sqrt((noncore_df['ca_dist_A']**2).mean()):.2f} Å  "
              f"max: {noncore_df['ca_dist_A'].max():.2f} Å")

    mobile = noncore_df[noncore_df['ca_dist_A'] > 5.0].sort_values('ca_dist_A', ascending=False)
    if not mobile.empty:
        print(f"Residues displaced > 5 Å (n={len(mobile)}): "
              f"{mobile['residue'].tolist()[:10]}{'...' if len(mobile) > 10 else ''}")
    return df

# Usage — non-repeat protein
df = analyze_conformational_change('apo.pdb', 'holo.pdb', q_threshold=3.0)
# df columns: residue, ca_dist_A, in_core
```

## Step 3: Workflow for Repeat Proteins

```python
# 1. Run symmetry analysis on both states to identify repeat units
#    (see symmetry-analysis.md — Internal Symmetry section)
#    Example: 7 repeat units of ~42 residues each

# 2. Extract corresponding repeat units from both states
def extract_domain(pdb_in, pdb_out, res_start, res_end):
    with open(pdb_in) as f, open(pdb_out, 'w') as g:
        for line in f:
            if line.startswith('ATOM') and res_start <= int(line[22:26]) <= res_end:
                g.write(line)
        g.write('END\n')

unit_size = 40  # from symmetry analysis — replace with your actual repeat unit size
for i in range(8):
    r1, r2 = 1 + i * unit_size, (i + 1) * unit_size
    extract_domain('state1.pdb', f'/tmp/s1_unit{i}.pdb', r1, r2)
    extract_domain('state2.pdb', f'/tmp/s2_unit{i}.pdb', r1, r2)

# 3. Analyze each unit independently
import pandas as pd
results = []
for i in range(8):
    df_unit = analyze_conformational_change(
        f'/tmp/s1_unit{i}.pdb', f'/tmp/s2_unit{i}.pdb',
        q_threshold=3.0
    )
    df_unit['unit'] = i
    results.append(df_unit)

df_all = pd.concat(results, ignore_index=True)
print(df_all.groupby('unit')['ca_dist_A'].agg(['mean', 'max']))
```

## Step 4: Interpreting Results

- **`in_core=True` residues** — structurally invariant anchor; displacement here reflects
  alignment noise, not true conformational change
- **`in_core=False` with large displacement** — regions that genuinely move between states
- **Large jump between core and non-core RMSD** (e.g., core 0.5 Å vs non-core 10+ Å)
  → rigid-body domain movement
- **Uniform small displacement everywhere** → subtle conformational change or no change

## Step 5 (Optional): Visualization

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
