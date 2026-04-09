# MICAN Symmetry & Domain Analysis

## Table of Contents
1. [Internal Symmetry — Basic Detection](#1-internal-symmetry--basic-detection)
2. [Rotation Matrix Analysis](#2-rotation-matrix-analysis)
3. [Cyclic Permutation — Algebraic Proof](#3-cyclic-permutation--algebraic-proof)
4. [Domain Detection from Sub-optimal Alignments](#4-domain-detection-from-sub-optimal-alignments)
5. [Distinguishing True Symmetry from False Positives](#5-distinguishing-true-symmetry-from-false-positives)

---

## 1. Internal Symmetry — Basic Detection

Compare a protein against itself. Rank 1 is always perfect self-match; Rank 2+ expose
internal symmetry or repeated domains.

> **⚠ Use `-g 1000` for self-comparison** — internal symmetry relies on sub-optimal solutions,
> so thorough search coverage is essential. Use `-n 20` to enumerate multiple solutions.

```python
import subprocess
from pymican import BINFILEPATH

pdb = 'protein.pdb'

result = subprocess.run(
    [BINFILEPATH, pdb, pdb, '-w', '-g', '1000', '-n', '20'],
    capture_output=True, text=True
)
in_table = False
for line in result.stdout.split('\n'):
    if 'Brief description' in line: in_table = True
    if in_table: print(line)
    if in_table and '(TMscore was' in line: break
# DALI Z > 2 in Rank 2+ → statistically significant internal symmetry
```

To get residue-level pairs for a specific rank:

```python
def get_self_alignment(pdb, rank, options='-w -g 1000 -n 20'):
    result = subprocess.run(
        [BINFILEPATH, '-z', pdb, pdb] + options.split() + ['-i', str(rank)],
        capture_output=True, text=True
    )
    rmsd = length = tmscore = dali = ''
    pairs = []
    for line in result.stdout.split('\n'):
        if line.startswith('# RMSD'):           rmsd = line.split()[2]
        elif line.startswith('# Length'):       length = line.split()[2]
        elif line.startswith('# TM-score: P1'): tmscore = line.split()[3]
        elif line.startswith('# Dali Zscore'):  dali = line.split()[2]
        elif line.startswith('ALIGN') and line.split()[4] != '.':
            pairs.append((int(line.split()[1]), int(line.split()[4])))
    return {'rmsd': rmsd, 'length': length, 'tmscore': tmscore, 'dali': dali, 'pairs': pairs}
```

---

## 2. Rotation Matrix Analysis

Extract rotation angle and axis from each sub-optimal alignment to identify the symmetry type.

Key equations:
- Rotation angle: `cos θ = (trace(R) − 1) / 2`
- Rotation axis: eigenvector of R with eigenvalue = 1
- Screw pitch `t∥ = t · axis` — if ≈ 0: pure rotation (Cn); if ≠ 0: screw symmetry

If multiple ranks share the same rotation axis (< 1° deviation), a true symmetry axis exists.
Rotation angles matching `k/n × 360°` (within ~1°) confirm Cn symmetry.

```python
import numpy as np
import subprocess
from pymican import BINFILEPATH

def get_rot_vec_pairs(pdb, rank, options='-w -g 1000 -n 20'):
    result = subprocess.run(
        [BINFILEPATH, '-z', pdb, pdb] + options.split() + ['-i', str(rank)],
        capture_output=True, text=True
    )
    rot = np.zeros((3, 3)); vec = np.zeros(3); pairs = []
    for line in result.stdout.split('\n'):
        if line.startswith(' 1   '):
            _, vec[0], rot[0,0], rot[0,1], rot[0,2] = map(float, line.split())
        if line.startswith(' 2   '):
            _, vec[1], rot[1,0], rot[1,1], rot[1,2] = map(float, line.split())
        if line.startswith(' 3   '):
            _, vec[2], rot[2,0], rot[2,1], rot[2,2] = map(float, line.split())
        if line.startswith('ALIGN') and line.split()[4] != '.':
            pairs.append((int(line.split()[1]), int(line.split()[4])))
    return rot, vec, pairs

def analyze_rotation(rot, vec):
    cos_theta = np.clip((np.trace(rot) - 1) / 2, -1, 1)
    theta_deg = np.degrees(np.arccos(cos_theta))
    eigvals, eigvecs = np.linalg.eig(rot)
    axis = eigvecs[:, np.argmin(np.abs(eigvals - 1.0))].real
    axis = axis / np.linalg.norm(axis)
    if axis[np.argmax(np.abs(axis))] < 0: axis = -axis
    t_par  = np.dot(vec, axis)
    t_perp = np.linalg.norm(vec - t_par * axis)
    return {'theta_deg': theta_deg, 'axis': axis, 't_par': t_par, 't_perp': t_perp,
            'det': np.linalg.det(rot)}

# Analyze Rank 2–8 and check axis alignment
pdb = 'protein.pdb'
results = {}
print(f"{'Rank':>5}  {'Angle':>8}  {'t∥(Å)':>7}  {'t⊥(Å)':>7}")
for rank in range(2, 9):
    rot, vec, pairs = get_rot_vec_pairs(pdb, rank)
    r = analyze_rotation(rot, vec)
    results[rank] = r
    print(f"{rank:>5}  {r['theta_deg']:>8.2f}  {r['t_par']:>7.3f}  {r['t_perp']:>7.2f}")

# Check axis alignment across ranks
for j in range(3, 9):
    dot = abs(np.dot(results[2]['axis'], results[j]['axis']))
    diff = np.degrees(np.arccos(np.clip(dot, 0, 1)))
    shared = '→ shared axis ✓' if diff < 1.0 else f'diff={diff:.1f}°'
    print(f"  Rank 2 vs {j}: axis diff = {diff:.2f}°  {shared}")

# Test Cn hypothesis
for n in [2, 3, 4, 6, 7, 8]:
    base = 360 / n
    matches = sum(1 for r in results.values()
                  if min(abs(r['theta_deg'] - base*k) for k in range(1, n)) < 2.0)
    print(f"  C{n}: {matches}/{len(results)} angles match (base={base:.1f}°)")

# Detect pseudo higher-order symmetry (e.g. pseudo C(2n) when true symmetry is Cn)
# Set TRUE_N to the confirmed symmetry order before running this block.
TRUE_N = 4  # <-- set to confirmed Cn order
PSEUDO_N = TRUE_N * 2
pseudo_base = 360 / PSEUDO_N
PSEUDO_THRESHOLD = 10.0  # degrees — loose threshold for pseudo-symmetry match

true_multiples = {360 / TRUE_N * k for k in range(1, TRUE_N)}
ref_axis = results[min(results)]['axis']  # reference axis from lowest-rank true-Cn result
pseudo_candidates = []

for rank, r in results.items():
    axis_diff = np.degrees(np.arccos(np.clip(abs(np.dot(r['axis'], ref_axis)), 0, 1)))
    if axis_diff > 1.0 or abs(r['t_par']) > 1.0:
        continue  # different axis or screw component — skip
    if any(abs(r['theta_deg'] - m) < 2.0 for m in true_multiples):
        continue  # already accounted for by true Cn
    dist = min(abs(r['theta_deg'] - pseudo_base * k) for k in range(1, PSEUDO_N))
    if dist < PSEUDO_THRESHOLD:
        pseudo_candidates.append((rank, r['theta_deg'], dist))

if pseudo_candidates:
    print(f"\n→ Pseudo C{PSEUDO_N} character detected (Δ < {PSEUDO_THRESHOLD}°):")
    for rank, angle, diff in pseudo_candidates:
        print(f"  Rank {rank}: {angle:.1f}°  (Δ{diff:.1f}° from exact C{PSEUDO_N} angle)")
    print(f"  Conclusion: True C{TRUE_N} symmetry.")
    print(f"  Visually pseudo C{PSEUDO_N} — repeat units are structurally similar but not")
    avg_dev = sum(d for _, _, d in pseudo_candidates) / len(pseudo_candidates)
    print(f"  strictly equivalent (average angle deviation ~{avg_dev:.1f}° from exact).")
else:
    print(f"\n→ No pseudo C{PSEUDO_N} character detected.")
```

---

## 3. Cyclic Permutation — Algebraic Proof

For tandem-repeat proteins, map each residue to its repeat unit and verify that unit-level
correspondence forms a cyclic permutation. This algebraically proves rotational symmetry
independent of geometric analysis.

- Cyclic permutation of units ↔ closed-ring arrangement ↔ Cn rotational symmetry
- All cyclic shifts {0, 1, ..., n−1} present across ranks → the full cyclic group Cn confirmed

The repeat unit size is estimated from the GCD of alignment shifts:

```python
from math import gcd
from functools import reduce
from collections import Counter

def estimate_unit_size(all_pairs_by_rank):
    """Estimate repeat unit length from alignment shifts.

    Strategy: compute GCD only from the lowest-rank sub-optimal alignments
    (Rank 2 and 3), which represent the dominant symmetry operations with the
    cleanest shifts. Higher ranks mix multiple symmetry operations and
    introduce fractional shifts that collapse the GCD to 1.

    The lowest rank with a non-trivial shift (≠ 0) gives the most reliable
    estimate of the fundamental repeat unit.
    """
    # Use only the first few ranks which have the cleanest shifts
    low_rank_shifts = []
    for rank in sorted(all_pairs_by_rank.keys())[:3]:  # Rank 2, 3, 4
        pairs = all_pairs_by_rank[rank]
        low_rank_shifts.extend(abs(r2 - r1) for r1, r2 in pairs if r2 != r1)

    if low_rank_shifts:
        result = reduce(gcd, low_rank_shifts)
        if result > 1:
            return result

    # Fallback: use the dominant (most common) non-zero shift across all ranks
    all_shifts = [abs(r2 - r1)
                  for pairs in all_pairs_by_rank.values()
                  for r1, r2 in pairs if r2 != r1]
    if not all_shifts:
        return None
    dominant = Counter(all_shifts).most_common(1)[0][0]
    return dominant

def check_cyclic_permutation(pairs, unit_size, start_res):
    """Return unit mapping and cyclic shift index (None if not cyclic)."""
    def get_unit(res): return (res - start_res) // unit_size
    mapping = {get_unit(r1): get_unit(r2) for r1, r2 in pairs}
    n = max(mapping.keys()) + 1
    ordered = [mapping.get(u) for u in range(n)]
    shifts = [[list(range(n))[(i + k) % n] for i in range(n)] for k in range(n)]
    shift_idx = next((k for k, s in enumerate(shifts) if ordered == s), None)
    return ordered, shift_idx, n

# Example workflow
pdb = 'protein.pdb'

# Step 1: get rotation results and pairs for all ranks
rotation_results = {}
all_pairs = {}
for rank in range(2, 9):
    rot, vec, pairs = get_rot_vec_pairs(pdb, rank)
    rotation_results[rank] = analyze_rotation(rot, vec)
    all_pairs[rank] = pairs

# Step 2: filter out false positives (rotation angle < 20°)
true_ranks = {rank for rank, r in rotation_results.items() if r['theta_deg'] > 30.0}
print(f"True symmetry ranks (angle > 30°): {sorted(true_ranks)}")

# Step 3: estimate unit size (uses lowest-rank shifts for robustness)
unit_size = estimate_unit_size(all_pairs)
start_res = min(r1 for pairs in all_pairs.values() for r1, _ in pairs)
print(f"Estimated unit_size={unit_size}  start_res={start_res}")

shift_set = set()
for rank, pairs in all_pairs.items():
    ordered, shift, n_units = check_cyclic_permutation(pairs, unit_size, start_res)
    shift_set.add(shift)
    print(f"Rank {rank}: shift={shift}  {[f'u{o}' for o in ordered]}")

print(f"Cyclic group C{n_units} confirmed: {shift_set == set(range(n_units))}")
```

---

## 4. Domain Detection from Sub-optimal Alignments

When a protein consists of multiple structural domains, sub-optimal alignments reveal
characteristic signals. **Look for these indicators before concluding Cn symmetry:**

### Signals that suggest a multi-domain protein

| Signal | Observation | Interpretation |
|--------|-------------|----------------|
| Rank 2–K: tight residue range | res1 max ≤ X for all high-scoring ranks | Repeat domain ends at ~X |
| Rank K+1+: sudden range expansion | res1 max jumps to full chain | Additional domain detected |
| Non-Cn rotation angles | Some ranks show angles not matching k/n × 360° | Those ranks span domain boundary |
| Mixed shift values | Shifts inconsistent with a single unit_size | Multiple structural units present |

### Workflow: detect domain boundary from sub-optimal alignments

```python
import subprocess
import numpy as np
from pymican import BINFILEPATH

pdb = 'protein.pdb'

# Step 1: collect residue range and rotation angles for each rank
print("Rank  r1_max  r2_max  theta  shifts")
all_data = {}
for rank in range(2, 16):
    result = subprocess.run(
        [BINFILEPATH, '-z', pdb, pdb, '-w', '-g', '1000', '-n', '20', '-i', str(rank)],
        capture_output=True, text=True
    )
    pairs = [(int(l.split()[1]), int(l.split()[4]))
             for l in result.stdout.split('\n')
             if l.startswith('ALIGN') and l.split()[4] != '.']
    rot = np.zeros((3,3))
    for l in result.stdout.split('\n'):
        if l.startswith(' 1   '): _, _, rot[0,0], rot[0,1], rot[0,2] = map(float, l.split())
        if l.startswith(' 2   '): _, _, rot[1,0], rot[1,1], rot[1,2] = map(float, l.split())
        if l.startswith(' 3   '): _, _, rot[2,0], rot[2,1], rot[2,2] = map(float, l.split())
    theta = np.degrees(np.arccos(np.clip((np.trace(rot)-1)/2, -1, 1)))
    shifts = sorted(set(r2-r1 for r1,r2 in pairs))[:4] if pairs else []
    r1_max = max((r1 for r1,_ in pairs), default=0)
    r2_max = max((r2 for _,r2 in pairs), default=0)
    all_data[rank] = {'pairs': pairs, 'theta': theta, 'r1_max': r1_max, 'r2_max': r2_max}
    print(f"{rank:>4}  {r1_max:>6}  {r2_max:>6}  {theta:>6.1f}°  {shifts}")

# Step 2: identify the boundary — where r1_max jumps
r1_maxes = [(rank, d['r1_max']) for rank, d in all_data.items()]
boundary_rank = next((rank for rank, mx in r1_maxes
                      if mx > r1_maxes[0][1] * 1.15), None)
if boundary_rank:
    prev_max = all_data[boundary_rank-1]['r1_max']
    curr_max = all_data[boundary_rank]['r1_max']
    print(f"\n→ Domain boundary detected near residue {prev_max}–{curr_max} "
          f"(Rank {boundary_rank-1}→{boundary_rank})")
```

### Analyzing each domain separately

Once a boundary is identified, extract each domain and analyze independently:

```python
def extract_domain(pdb_in, pdb_out, res_start, res_end):
    with open(pdb_in) as f, open(pdb_out, 'w') as g:
        for line in f:
            if line.startswith('ATOM') and res_start <= int(line[22:26]) <= res_end:
                g.write(line)
        g.write('END\n')

# Split and analyze each domain
boundary = 295  # example boundary
extract_domain(pdb, '/tmp/domain1.pdb', 1, boundary)
extract_domain(pdb, '/tmp/domain2.pdb', boundary + 1, 999)

# Self-comparison of each domain
for name, dp in [('Domain1', '/tmp/domain1.pdb'), ('Domain2', '/tmp/domain2.pdb')]:
    result = subprocess.run(
        [BINFILEPATH, dp, dp, '-w', '-g', '1000', '-n', '5'],
        capture_output=True, text=True
    )
    print(f"\n{name} self-comparison:")
    in_table = False
    for line in result.stdout.split('\n'):
        if 'Brief description' in line: in_table = True
        if in_table: print(f"  {line}")
        if in_table and '(TMscore was' in line: break

# Inter-domain comparison
from pymican import mican
m = mican()
aln = m.align('/tmp/domain1.pdb', '/tmp/domain2.pdb', options='-w')
tm_mean = (aln.TMscore1 + aln.TMscore2) / 2
print(f"\nDomain1 vs Domain2: TMscore={aln.TMscore:.3f}  TM_mean={tm_mean:.3f}  "
      f"RMSD={aln.rmsd:.2f}Å  nalign={aln.nalign}  SeqID={aln.seq_identity:.1f}%")
print(f"  coverage1={aln.coverage1:.1f}%  coverage2={aln.coverage2:.1f}%")
```

---

## 5. Distinguishing True Symmetry from False Positives

Not all sub-optimal solutions with DALI Z > 2 represent true symmetry operations.

### True symmetry operations

- Rotation angle > 30°, matching k/n × 360° within ~1°
- Correspondence between **different** repeat units (u0→u1, u0→u2, etc.)
- Cyclic permutation pattern across multiple ranks
- Shared rotation axis (< 1° deviation) across all matching ranks
- t∥ ≈ 0 (pure rotation, no screw component)

### False positives from spatial proximity

- Very small rotation angle (< 20°) — nearly identity transformation
- Correspondence **within the same unit** (u0→u0)
- No cyclic permutation
- Often appear when a C-terminal domain is spatially close to the repeat region

> A rank with rotation angle < 20° and same-unit correspondence should be excluded from
> symmetry analysis. These reflect packing geometry, not molecular symmetry.
