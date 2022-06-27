from email.policy import default
import os
import tqdm
import numpy as np
import pandas as pd
from collections import defaultdict
from pymican import mican

# Dataset
dataset_name = 'MALISAM-ns'
dir_tests = os.path.dirname(os.path.realpath(__file__))
listfile = dir_tests+'/dataset/'+dataset_name+'/list'
with open(listfile, 'r') as f:
    targets = [l.strip()  for l in f.readlines()]

# MICAN object
m = mican()

stack = []
for t in tqdm.tqdm(targets):
    pdbfile1 = dir_tests+'/dataset/'+dataset_name+'/pdb/'+t+'.1.pdb'
    pdbfile2 = dir_tests+'/dataset/'+dataset_name+'/pdb/'+t+'.2.pdb'
    # alignment
    aln_out = m.align(pdbfile1, pdbfile2, options='-w')['alignment']
    # alignment identity check
    aln_ref = dir_tests+'/dataset/'+dataset_name+'/csv/'+t+'.csv'
    aln_ref = pd.read_csv(aln_ref, index_col=0)
    assert all(
        np.array(aln_out['residue1'].values, dtype=int)
        == np.array(aln_ref['residue1'].values, dtype=int)
        ), print("Unequal resutls in {}.".format(t))
    # to dict    
    aln_dict = defaultdict(lambda: None)
    for i,j in zip(aln_out['residue1'].values, aln_out['residue2'].values):
        aln_dict[int(i)] = int(j)
    # reference alignment
    alnfile = dir_tests+'/dataset/'+dataset_name+'/reference/'+t+'.aln'
    with open(alnfile, 'r') as f:
        aln_targets = np.array([list(map(int, l.split()))  for l in f.readlines()])
    # check accuracy
    check = np.array([aln_dict[a[0]]==a[1] for a in aln_targets])
    stack.append(check.mean())

accuracy_mean = np.array(stack).mean()
assert accuracy_mean >= 0.676, 'Low accuracy in [{}]'.format(dataset_name)
print("accuracy_mean = {:.3f} >= 0.676.".format(accuracy_mean))