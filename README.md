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
MIT

# Easy instllation
```
pip install pymican
```

# Compilation and usage
1. To compile MICAN software: please type following command
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

For more details, please read following usage.

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

# Usage as Python module
1. install pymican
```
pip install pymican
```
2. usage
```python
from pymican import mican

m = mican()

outdict = m.align(pdb1=INPUT_PDBFILE_1, pdb2=IPUT_PDBFILE_2, options=EXTRAOPTIONS)
```

```python
Paremeters
  ----------
  pdb1, pdb2 : str
      Input PDB files
  options : (str, [str,...]), default=[]
      Extra potions for mican calculation.
      For the option details please see (https://github.com/ShintaroMinami/mican).
  Returns
  -------
  dict = {
      'mode': ('sequential', 'rewirering', 'reverse'),
      'pdb1': string,
      'size1': int,
      'pdb2': string,
      'size2': int,
      'nalign': int,
      'rmsd': float,
      'TMscore': float,
      'sTMscore': float,
      'seq-identity': float,
      'DALIscore': float,
      'SPscore': float,
      'TMscore1': float, # TMscore normalized by size of protein1
      'coverage1': float,
      'TMscore2': float, # TMscore normalized by size of protein2
      'coverage2': float
  }
```
