import re

def output2dict(mican_output:str):
    dict = {}
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
            dict['seq-identity'] = float(seqid)
        if l.startswith(' Alignment mode'):
            dict['mode'] = l.split()[3]
        if l.startswith(' Protein1'):
            dict['size1'] = int(re.findall('(....) residues', l)[0])
            dict['pdb1'] = str(l.split('= ')[1])
        if l.startswith(' Protein2'):
            dict['size2'] = int(re.findall('(....) residues', l)[0])
            dict['pdb2'] = str(l.split('= ')[1])
    return dict