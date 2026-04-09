"""
MICAN alignment script to find residue correspondence.
Aligns two PDB files using non-sequential mode and extracts correspondence for residues 10, 20, 30.
"""
from pymican import mican
import pandas as pd
import os

# Define PDB files
pdb1 = "/Users/sminami/github/mican/tests/dataset/MALISAM-ns/pdb/d1a2za_d1ghha_.1.pdb"
pdb2 = "/Users/sminami/github/mican/tests/dataset/MALISAM-ns/pdb/d1a2za_d1ghha_.2.pdb"

# Output directory
output_dir = "/Users/sminami/github/mican/.claude/skills/mican-analysis-workspace/iteration-1/eval-2/with_skill/outputs"

# Initialize MICAN
m = mican()

# Perform non-sequential alignment (-w option)
output_lines = []
output_lines.append("Running MICAN alignment with non-sequential mode (-w)...")
output_lines.append(f"PDB1: {pdb1}")
output_lines.append(f"PDB2: {pdb2}")
output_lines.append("")

aln = m.align(pdb1=pdb1, pdb2=pdb2, options="-w")

# Display alignment scores
output_lines.append("=" * 60)
output_lines.append("Alignment Results")
output_lines.append("=" * 60)
output_lines.append(f"TM-score:         {aln.TMscore:.4f}")
output_lines.append(f"RMSD:             {aln.rmsd:.2f} Å")
output_lines.append(f"Aligned residues: {aln.nalign}/{aln.size1}")
output_lines.append(f"Sequence identity: {aln.seq_identity:.1f}%")
output_lines.append("")

# Target residues to query
target_residues = [10, 20, 30]

# Extract correspondence information
df = aln.alignment.copy()
df["residue1"] = df["residue1"].astype(str)

output_lines.append("=" * 60)
output_lines.append("Residue Correspondence Table")
output_lines.append("=" * 60)

rows = []
for res in target_residues:
    hit = df[df["residue1"] == str(res)]
    if not hit.empty:
        r = hit.iloc[0]
        rows.append({
            "pdb1_residue": res,
            "pdb1_chain": r["chain1"],
            "pdb1_aa": r["aatype1"],
            "pdb2_residue": r["residue2"],
            "pdb2_chain": r["chain2"],
            "pdb2_aa": r["aatype2"],
            "ca_distance_A": float(r["distance"])
        })
        output_lines.append(f"PDB1 residue {res} ({r['chain1']}:{r['aatype1']}) → PDB2 residue {r['residue2']} ({r['chain2']}:{r['aatype2']}), Cα distance = {float(r['distance']):.2f} Å")
    else:
        rows.append({
            "pdb1_residue": res,
            "pdb1_chain": "-",
            "pdb1_aa": "-",
            "pdb2_residue": "N/A",
            "pdb2_chain": "-",
            "pdb2_aa": "-",
            "ca_distance_A": None
        })
        output_lines.append(f"PDB1 residue {res}: Not aligned (no correspondence)")

output_lines.append("")

# Create DataFrame
result_df = pd.DataFrame(rows)

# Save to CSV
output_csv = os.path.join(output_dir, "residue_correspondence.csv")
result_df.to_csv(output_csv, index=False)
output_lines.append(f"Results saved to: {output_csv}")
output_lines.append("")
output_lines.append("CSV content:")
output_lines.append(result_df.to_string(index=False))

# Print all output
full_output = "
".join(output_lines)
print(full_output)

# Save to result.txt
result_path = os.path.join(output_dir, "result.txt")
with open(result_path, "w") as f:
    f.write(full_output)
