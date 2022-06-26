# MICAN
Protein structure alignment program that can handle
- M: multiple-chain complexs
- I: inverse direction of SSEs
- C: Ca only models
- A: alternative alignments
- N: non-sequential alignments

## Author information
Author: S. Minami, K. Sawada, and G. Chikenji

Web Site: http://www.tbp.cse.nagoya-u.ac.jp/MICAN

## References
1. BMC Bioinformatics 2013, 14(24), S. Minami, K. Sawada, and G. Chikenji
2. Bioinformatics 2018, 34(19), S. Minami, K. Sawada, M Ota, and G. Chikenji

## License
[MIT](https://choosealicense.com/licenses/mit/)
# Easy instllation
```
pip install pymican
```

# Python module usage
1. install pymican
```
pip install pymican
```
2. usage
```python
# import module
from pymican import mican

# create object
m = mican()

# calculate alignment
aln = m.align(pdb1='pdbfile1', pdb2='pdbfile2', options='extra-mican-options')

# get TM-score, RMSD, etc.
print(aln.TMscore)
print(aln.rmsd)
```

Attributes of Alignment object
```
    MICAN alignment class

    Attributes
    ----------
    outdict : dict
        Alignment info
    mode : str
        Alignment mode
    pdb1, pdb2 : str
        PDB file path
    size1, size2 : int
        Size of protein structure
    nalign : int
        Number of aligned residues
    rmsd : float
        RMSD of aligned residues
    TMscore : float
        TM-score
    sTMscore : float
        SSE weighted TM-score
    seq_identity : float
        Sequence identity as percentage [0,100]
    DALIscore : float
        DALI z-score
    SPscore : float
        SP-score
    TMscore1, TMscore2 : float
        TM-score normalized by each protein length
    coverage1, coverage2 : float
        Aligned coverage for each protein length
    translation_rot : numpy.array(3,3)
        Rotation matrix for superposition protein1 on protein2
    translation_vec : numpy.array(3)
        Translation vector for superposition protein1 on protein2
    alignment : pandas.DataFrame
        Residue-Residue alignment info
    alignlst : List[pandas.item]
        Alignment info for iterator methods

    Methods
    -------
    translate_xyz(xyz: np.array(N,3)) -> np.array(N,3)
        Rotate & translate xyz coordinates
```

# Compilation and usage
1. To compile MICAN software: please type this command
```
% make
```

2. To run MICAN software:
```
% mican protein1 protein2 -a align.aln -o sup.pdb
```

--  e.g. --
```
% mican test/test1.1.pdb test/test1.2.pdb -a align.aln -o sup.pdb
```

For more details, please read the usage.

```
 USAGE: % mican protein1 protein2 [OPTION]

 Description:
  -f             fast mode (same as "-g 15")
  -s             sequential (SQ) alignment mode
  -w             rewiring (RW) alignment mode
  -r             rewiring & reverse (RR) alignment mode
  -R             reverse constrained alignment mode
  -x             silent mode (without any output on the console)
  -p             print alignment progress
  -c1 ChainIDs   chain ID specifier for protein1 (e.g. -c1 A, -c1 ABC)
  -c2 ChainIDs   chain ID specifier for protein2
  -o  Filename   superposition file (rasmol-script)
  -a  Filename   alignment file
  -m  Filename   translation matrix file
  -n  Integer    number of solutions output (default=5)
  -i  Integer    output i-th solution on stdout & superposition file
  -t  Integer    selection score ([0]:sTMscore, 1:TMscore, 2:SPscore)
  -g  Integer    number of GH candidates used (default=50)
  -l  Integer    minimum segment length (default=3)
  -d  Real       fix TM-score scaling factor d0
  -q  Real       maximum distance between Ca atoms to be aligned (default=10.0)
  
 Simple usage (SQ):
   % mican protein1 protein2
   % mican protein1 protein2 -a align.aln -o sup.pdb

 Rewiring mode alignment (RW):
   % mican protein1 protein2 -w

 Rewiring & reverse mode alignment (RR):
   % mican protein1 protein2 -r

 To visualize superposition:
   % mican protein1 protein2 -o sup.pdb
   % rasmol -script sup.pdb
```

