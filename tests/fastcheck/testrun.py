#! /usr/bin/env python

import sys
sys.path.append("../../")
import numpy as np
# import mican class
from pymican import mican

# setup mican object
mican = mican()

# calculate mican alignment for two pdb files
aln = mican.align('test1.1.pdb', 'test1.2.pdb')

# calculate mican alignment for pdb line and pdb file
pdbline1 = "".join(open('test1.1.pdb', 'r').readlines())
pdbline2 = "".join(open('test1.2.pdb', 'r').readlines())
aln = mican.align(pdbline1, pdbline2)

# print summary
print(aln)

# print TM-score
print(aln.TMscore)

# print the first 10 aligned residue-pairs
for a in aln[0:10]:
    print(a.residue1, a.residue2, a.distance)

# translate coordinates by superposition matrix
xyz = np.random.rand(10,3)
print(xyz)
xyz_translated = aln.translate_xyz(xyz)
print(xyz_translated)