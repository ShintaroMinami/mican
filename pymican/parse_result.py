import re
import numpy as np
import pandas as pd

def output2dict(mican_output:str):
    dict = {}
    rot = np.zeros((3,3), dtype=float)
    vec = np.zeros((3), dtype=float)
    align = []
    for l in mican_output.split('\n'):
        if l.startswith(' TM-score=') & l.endswith('Protein1)'):
            dict['TMscore1'], dict['coverage1'] = map(float, re.findall('=(.....)', l))
        if l.startswith(' TM-score=') & l.endswith('Protein2)'):
            dict['TMscore2'], dict['coverage2'] = map(float, re.findall('=(.....)', l))
        if l.startswith('    1    '):
            _, stm, tm, dali, sp, naln, rmsd, seqid = l.split()
            dict['sTMscore'] = float(stm)
            dict['TMscore'] = float(tm)
            dict['DALIscore'] = float(dali)
            dict['SPscore'] = float(sp)
            dict['nalign'] = int(naln)
            dict['rmsd'] = float(rmsd)
            dict['seq_identity'] = float(seqid)
        if l.startswith(' Alignment mode'):
            dict['mode'] = l.split()[3]
        if l.startswith(' Protein1'):
            dict['size1'] = int(re.findall('(....) residues', l)[0])
            dict['pdb1'] = str(l.split('= ')[1])
        if l.startswith(' Protein2'):
            dict['size2'] = int(re.findall('(....) residues', l)[0])
            dict['pdb2'] = str(l.split('= ')[1])
        if l.startswith(' 1   '):
            _, vec[0], rot[0,0], rot[0,1], rot[0,2] = [float(v) for v in l.split()]
        if l.startswith(' 2   '):
            _, vec[1], rot[1,0], rot[1,1], rot[1,2] = [float(v) for v in l.split()]
        if l.startswith(' 3   '):
            _, vec[2], rot[2,0], rot[2,1], rot[2,2] = [float(v) for v in l.split()]
        if l.startswith('ALIGN'):
            _, iaa1, aa1, ch1, iaa2, aa2, ch2, close, dist = l.split()
            if iaa2 != '.':
                align.append((iaa1, iaa2, ch1, ch2, aa1, aa2, close, dist))
    dict['translation_rot'] = rot
    dict['translation_vec'] = vec
    dict['alignment'] = pd.DataFrame(align,
        columns=['residue1', 'residue2','chain1', 'chain2',
                'aatype1', 'aatype2', 'closeness', 'distance'])
    # return
    return dict