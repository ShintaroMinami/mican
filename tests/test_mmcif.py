"""Check that PDBx/mmCIF input gives the same results as PDB input.

PDB files of the MALISAM-ns dataset are converted to mmCIF on the fly, and
the MICAN results for both formats must be identical.
Run directly (python tests/test_mmcif.py) or with pytest.
"""
import os
import tempfile

from pymican import mican

DIR_TESTS = os.path.dirname(os.path.realpath(__file__))
DATASET = os.path.join(DIR_TESTS, "dataset", "MALISAM-ns")
N_PAIRS = 10
ITEMS = [
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
    "label_comp_id", "label_asym_id", "label_seq_id", "pdbx_PDB_ins_code",
    "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv",
    "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id",
    "pdbx_PDB_model_num",
]


def _cif_value(v: str) -> str:
    if not v:
        return "?"
    return f'"{v}"' if "'" in v else v


def pdb_to_cif(pdb_text: str) -> str:
    """Convert the ATOM/HETATM records of a PDB file to mmCIF.

    Args:
        pdb_text: Contents of a PDB file.

    Returns:
        mmCIF text with an _atom_site loop (label_seq_id is left unknown).
    """
    rows, model = [], 1
    for line in pdb_text.splitlines():
        if line.startswith("MODEL"):
            model = int(line[10:14])
        if not line.startswith(("ATOM ", "HETATM")):
            continue
        line = line.ljust(80)
        name = line[12:16]
        element = line[76:78].strip() or (
            name[1] if name[0] in " 0123456789" else name[:2].strip())
        chain = line[21].strip() or "."
        resname = line[17:20].strip()
        rows.append([
            line[0:6].strip(), line[6:11].strip(), element, name.strip(),
            line[16].strip() or ".", resname, chain, ".",
            line[26].strip() or "?", line[30:38].strip(), line[38:46].strip(),
            line[46:54].strip(), line[54:60].strip(), line[60:66].strip(),
            line[22:26].strip(), resname, chain, name.strip(), str(model),
        ])
    header = "data_test\n#\nloop_\n" + "".join(f"_atom_site.{i}\n" for i in ITEMS)
    body = "".join(" ".join(_cif_value(v) for v in r) + "\n" for r in rows)
    return header + body + "#\n"


def test_mmcif_matches_pdb() -> None:
    """MICAN gives identical alignments and scores for PDB and mmCIF input."""
    m = mican()
    with open(os.path.join(DATASET, "list")) as f:
        targets = [line.strip() for line in f][:N_PAIRS]
    with tempfile.TemporaryDirectory() as tmp:
        for t in targets:
            pdbs, cifs = [], []
            for k in ("1", "2"):
                pdb = os.path.join(DATASET, "pdb", f"{t}.{k}.pdb")
                cif = os.path.join(tmp, f"{t}.{k}.cif")
                with open(pdb) as fin, open(cif, "w") as fout:
                    fout.write(pdb_to_cif(fin.read()))
                pdbs.append(pdb)
                cifs.append(cif)
            for opt in ("-w", "-s"):
                a = m.align(*pdbs, options=opt)
                b = m.align(*cifs, options=opt)
                assert a.alignment.equals(b.alignment), f"{t} {opt}: alignment"
                assert (a.TMscore, a.rmsd, a.DALIscore) == (
                    b.TMscore, b.rmsd, b.DALIscore), f"{t} {opt}: scores"
    print(f"mmCIF input matches PDB input ({len(targets)} pairs x 2 modes).")


if __name__ == "__main__":
    test_mmcif_matches_pdb()
