"""Helpers for MICAN results that pymican does not expose.

pymican's `mican().align()` returns only the first-ranked solution, hides
MICAN's error messages, and does not parse the per-chain results. These
helpers run the MICAN binary installed with pymican and parse its output.

Usage:
    import sys
    sys.dont_write_bytecode = True   # keep the skill directory free of __pycache__
    sys.path.insert(0, ".claude/skills/mican-analysis/scripts")
    from mican_utils import run_mican, brief_table, chain_summary, align_rank, write_superposed

    out = run_mican("a.pdb", "b.pdb", "-w -n 10 -g 500")
    print(brief_table(out))        # all ranks with their scores
    print(chain_summary(out))      # C3score and per-chain coverage (rank shown)
    aln = align_rank("a.pdb", "b.pdb", rank=2, options="-w -n 10 -g 500")
    write_superposed("a.pdb", "b.pdb", "a_on_b.pdb", options="-w")
"""
import os
import re
import subprocess
import tempfile
from typing import Sequence, Union

import pandas as pd
from pymican import BINFILEPATH
from pymican.main import Alignment
from pymican.parse_result import output2dict

Options = Union[str, Sequence[str]]


def _split(options: Options) -> list[str]:
    return options.split() if isinstance(options, str) else list(options)


def run_mican(pdb1: str, pdb2: str, options: Options = "-w",
              binary: str = BINFILEPATH) -> str:
    """Run MICAN and return its standard output.

    Args:
        pdb1: First structure file (PDB or mmCIF).
        pdb2: Second structure file (PDB or mmCIF).
        options: MICAN options, e.g. "-w -n 10".
        binary: MICAN executable (default: the one bundled with pymican).

    Returns:
        MICAN's standard output.

    Raises:
        RuntimeError: MICAN failed; the message is MICAN's own error text.
    """
    r = subprocess.run([binary, pdb1, pdb2, *_split(options)],
                       capture_output=True, text=True)
    errors = [l for l in (r.stdout + r.stderr).splitlines() if l.startswith("!! Error")]
    if r.returncode != 0 or errors:
        raise RuntimeError("\n".join(errors) or r.stderr.strip()
                           or f"mican exited with status {r.returncode}")
    return r.stdout


def brief_table(output: str) -> pd.DataFrame:
    """Parse the 'Brief description' table (one row per ranked solution).

    Args:
        output: MICAN standard output.

    Returns:
        DataFrame with columns rank, sTMscore, TMscore, DALIscore (Dali Z),
        SPscore, nalign, rmsd, seq_identity. TMscore is normalized by
        Protein1, as in MICAN's table.
    """
    rows, inside = [], False
    for line in output.splitlines():
        if "Brief description" in line:
            inside = True
        elif inside and line.startswith(" (TMscore was"):
            break
        elif inside and re.match(r"^\s+\d+\s+-?\d", line):
            v = line.split()
            rows.append([int(v[0]), *map(float, v[1:5]), int(v[5]), float(v[6]), float(v[7])])
    return pd.DataFrame(rows, columns=["rank", "sTMscore", "TMscore", "DALIscore",
                                       "SPscore", "nalign", "rmsd", "seq_identity"])


def chain_summary(output: str) -> dict:
    """Parse C3score and the per-chain results of the detailed solution.

    Args:
        output: MICAN standard output (details are for rank 1, or for the
            rank given with -i).

    Returns:
        {"C3score": float, "chains": DataFrame} where each row is a chain
        of Protein1 with coverage > 0: chain1, size1, coverage1 (% of the
        chain aligned) and partners ({chain2: % of chain1 aligned to it}).
    """
    m = re.search(r"C3score=\s*([\d.]+)", output)
    rows = []
    for line in output.splitlines():
        hit = re.match(r"^ \[P1:(.) \(size\s*(\d+)\)\] Coverage=\s*([\d.]+)% \((.*)\)", line)
        if hit:
            partners = {c: float(p) for c, p in re.findall(r"P2:(.)\s+([\d.]+)%", hit.group(4))}
            rows.append([hit.group(1), int(hit.group(2)), float(hit.group(3)), partners])
    return {"C3score": float(m.group(1)) if m else None,
            "chains": pd.DataFrame(rows, columns=["chain1", "size1", "coverage1", "partners"])}


def align_rank(pdb1: str, pdb2: str, rank: int = 1, options: Options = "-w",
               binary: str = BINFILEPATH) -> Alignment:
    """pymican Alignment for the rank-th solution, with that rank's scores.

    With `-i N`, pymican's `align()` mixes ranks: alignment, TMscore1/2,
    coverage and the rotation are rank N, but TMscore, sTMscore, rmsd,
    nalign, DALIscore, SPscore and seq_identity are rank 1. This function
    takes those from rank N's row of the brief table instead.

    Args:
        pdb1: First structure file.
        pdb2: Second structure file.
        rank: Solution rank (1 = best).
        options: MICAN options; -n is raised to `rank` if needed.
        binary: MICAN executable.

    Returns:
        Alignment object for the requested rank.

    Raises:
        ValueError: MICAN found fewer than `rank` solutions (try a larger -g).
    """
    opts = _split(options)
    if "-n" in opts:
        i = opts.index("-n")
        opts[i + 1] = str(max(int(opts[i + 1]), rank))
    else:
        opts += ["-n", str(max(5, rank))]
    out = run_mican(pdb1, pdb2, ["-z", *opts, "-i", str(rank)], binary)
    table = brief_table(out).set_index("rank")
    if rank not in table.index:
        raise ValueError(f"MICAN found only {len(table)} solutions; increase -g (e.g. -g 500)")
    row = table.loc[rank]
    d = output2dict(out)
    d.update(sTMscore=row.sTMscore, TMscore=row.TMscore, DALIscore=row.DALIscore,
             SPscore=row.SPscore, nalign=int(row.nalign), rmsd=row.rmsd,
             seq_identity=row.seq_identity)
    return Alignment(d)


def write_superposed(pdb1: str, pdb2: str, out_pdb: str, options: Options = "-w",
                     binary: str = BINFILEPATH) -> str:
    """Write protein1 superposed onto protein2 as a plain PDB file.

    MICAN's -o file is a RasMol script (MODEL 1 = protein1 moved onto
    protein2, MODEL 2 = protein2 unchanged). This keeps MODEL 1 only, so the
    result can be opened together with the original pdb2 in any viewer.

    Args:
        pdb1: Structure to move.
        pdb2: Reference structure (its coordinates are not changed).
        out_pdb: Output PDB path.
        options: MICAN options (add "-i N" for the N-th solution).
        binary: MICAN executable.

    Returns:
        out_pdb.
    """
    with tempfile.TemporaryDirectory() as tmp:
        sup = os.path.join(tmp, "sup.ras")
        run_mican(pdb1, pdb2, [*_split(options), "-o", sup], binary)
        lines, model = [], None
        for line in open(sup):
            if line.startswith("MODEL"):
                model = int(line.split()[1])
            elif line.startswith("ENDMDL"):
                model = None
            elif model == 1 and line.startswith(("ATOM", "HETATM")):
                lines.append(line)
    with open(out_pdb, "w") as f:
        f.writelines(lines)
        f.write("END\n")
    return out_pdb
