#!/usr/bin/env python3
"""
Compare two PDB structures using MICAN RW mode and extract residue correspondence.
"""

import os
import pandas as pd
from pymican import mican

# Constants
CA_DIST_WARN = 4.5  # Threshold for structural correspondence reliability
PDB1 = "tests/dataset/MALISAM-ns/pdb/d1a2za_d1ghha_.1.pdb"
PDB2 = "tests/dataset/MALISAM-ns/pdb/d1a2za_d1ghha_.2.pdb"
TARGET_RESIDUES = [10, 20, 30]
OUTPUT_DIR = ".claude/skills/mican-analysis-workspace/iteration-3/eval-2/with_skill/outputs"

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Run MICAN alignment
print("Running MICAN alignment in RW mode...")
m = mican()
aln = m.align(pdb1=PDB1, pdb2=PDB2, options="-w")

# Print summary scores
print("
=== Alignment Summary ===")
print(f"TM-score: {aln.TMscore:.4f}")
print(f"RMSD:     {aln.rmsd:.2f} Å")
print(f"Aligned:  {aln.nalign}/{aln.size1} residues")
print(f"Seq ID:   {aln.seq_identity:.1f}%")

# Extract correspondence for target residues
print(f"
=== Residue Correspondence (pdb1 → pdb2) ===")
df = aln.alignment.copy()
df["residue1"] = df["residue1"].astype(str)

rows = []
for res in TARGET_RESIDUES:
    hit = df[df["residue1"] == str(res)]
    if not hit.empty:
        r = hit.iloc[0]
        dist = float(r["distance"])
        rows.append({
            "pdb1_residue": res,
            "pdb1_chain": r["chain1"],
            "pdb1_aa": r["aatype1"],
            "pdb2_residue": r["residue2"],
            "pdb2_chain": r["chain2"],
            "pdb2_aa": r["aatype2"],
            "ca_dist_A": round(dist, 2),
            "note": "⚠ 対応不確か (>4.5Å)" if dist > CA_DIST_WARN else "OK"
        })
    else:
        rows.append({
            "pdb1_residue": res,
            "pdb1_chain": "-",
            "pdb1_aa": "?",
            "pdb2_residue": "N/A",
            "pdb2_chain": "-",
            "pdb2_aa": "-",
            "ca_dist_A": None,
            "note": "アライメント領域外"
        })

# Create result table
table = pd.DataFrame(rows)

# Print to console
print(table.to_string(index=False))

# Save as CSV
csv_path = f"{OUTPUT_DIR}/residue_correspondence.csv"
table.to_csv(csv_path, index=False)
print(f"
CSV saved: {csv_path}")

# Check for warnings
warned = table[table["note"].str.startswith("⚠", na=False)]
if not warned.empty:
    print(f"
⚠ 注意: {len(warned)}残基でCα距離が{CA_DIST_WARN}Åを超えています。")
    print("これらの対応は構造的に離れており、機能的な対応とは言えない可能性があります。")

# Save detailed result
result_path = f"{OUTPUT_DIR}/result.txt"
with open(result_path, "w") as f:
    f.write("=== MICAN Alignment Result ===
")
    f.write(f"PDB1: {PDB1}
")
    f.write(f"PDB2: {PDB2}
")
    f.write(f"Mode: RW (汎用アライメント)

")
    f.write(f"TM-score: {aln.TMscore:.4f}
")
    f.write(f"RMSD:     {aln.rmsd:.2f} Å
")
    f.write(f"Aligned:  {aln.nalign}/{aln.size1} residues
")
    f.write(f"Seq ID:   {aln.seq_identity:.1f}%

")
    f.write("=== Residue Correspondence ===
")
    f.write(table.to_string(index=False))
    if not warned.empty:
        f.write(f"

⚠ Warning: {len(warned)} residue(s) exceed {CA_DIST_WARN}Å threshold.
")

print(f"Result saved: {result_path}")
