import os
import sys
import subprocess
from typing import List, Union
import numpy as np
import pandas as pd
from .parse_result import output2dict

if sys.platform == "linux":
    try:
        from memory_tempfile import MemoryTempfile
        tempfile = MemoryTempfile()
    except ImportError:
        import tempfile
else:
    import tempfile


dir_script = os.path.dirname(os.path.realpath(__file__))
BINFILEPATH = os.path.abspath(dir_script+'/bin/mican')


class Alignment:
    """
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
    """

    def __init__(self, outdict: dict):
        self.outdict = outdict
        self.mode = outdict['mode']
        self.pdb1 = outdict['pdb1']
        self.size1 = outdict['size1']
        self.pdb2 = outdict['pdb2']
        self.size2 = outdict['size2']
        self.nalign = outdict['nalign']
        self.rmsd = outdict['rmsd']
        self.TMscore = outdict['TMscore']
        self.sTMscore = outdict['sTMscore']
        self.seq_identity = outdict['seq_identity']
        self.DALIscore = outdict['DALIscore']
        self.SPscore = outdict['SPscore']
        self.TMscore1 = outdict['TMscore1']
        self.coverage1 = outdict['coverage1']
        self.TMscore2 = outdict['TMscore2']
        self.coverage2 = outdict['coverage2']
        self.translation_rot = outdict['translation_rot']
        self.translation_vec = outdict['translation_vec']
        self.alignment = outdict['alignment']
        self.alignlist = list(self.alignment.itertuples())

    def __iter__(self):
        return iter(self.alignlist)

    def __str__(self):
        return '{}'.format(self.outdict)

    def __getitem__(self, key):
        if type(key) == str:
            return self.outdict.get(key)
        elif type(key) == int:
            return self.alignlist[key]
        elif type(key) == slice:
            return self.alignlist[key]

    def keys(self):
        return self.outdict.keys()

    def translate_xyz(self, xyz: Union[np.ndarray, List[float]])->np.ndarray:
        """
        Translate coordinates

        Paremeters
        ----------
        xyz : (np.array(3) | np.array(N,3) | List[float])
            Input xyz coordinates
        
        Returns
        -------
        np.array
            Translated xyz coordinates
        """
        # list to np.array
        xyz = np.array(xyz, dtype=float) if type(xyz) == list else xyz
        # save original shape
        original_shape = xyz.shape
        # add dimention 1 if ndim==1
        xyz = xyz[np.newaxis,:] if xyz.ndim == 1 else xyz
        # translate
        xyz_translated = np.dot(xyz,self.translation_rot.T) + self.translation_vec
        # return
        return np.reshape(xyz_translated, original_shape)


class mican:
    """
    MICAN: non-sequential alignment algorithm

    Attributes
    ----------
    binary : Filename
        executable binary file path
    """
    def __init__(self, binary: str=BINFILEPATH):
        """
        Parameters
        ----------
        binary : Filename
            executable binary file path
        """
        self.binary = binary
        return

    def align(self, pdb1: str, pdb2: str, options: Union[str, List[str]]=[], return_dict: int=False)->Union[Alignment, dict]:
        """
        Alignment calculation

        Paremeters
        ----------
        pdb1, pdb2 : str
            Input PDB files.
        options : (str | [str,...]), default=[]
            Extra potions for mican calculation.
            For the option details please see (https://github.com/ShintaroMinami/mican).
        return_dict : bool, defualt=False
            Simple dict return.

        Returns
        -------
        Alignment object

        if (return_dict == True)
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
            'seq_identity': float,
            'DALIscore': float,
            'SPscore': float,
            'TMscore1': float, # TMscore normalized by size of protein1
            'coverage1': float,
            'TMscore2': float, # TMscore normalized by size of protein2
            'coverage2': float,
            'translation_rot': numpy.array(3,3),
            'translation_vec': numpy.array(3),
            'alignment': pandas.DataFrame
        }
        """
        # arguments
        options = options if type(options) == str else ' '.join(options)
        # pdb files
        ntf1, ntf2 = None, None
        if not os.path.isfile(pdb1):
            ntf1 = tempfile.NamedTemporaryFile(mode='w')
            ntf1.write(pdb1)
            ntf1.flush()
            pdb1 = ntf1.name
        if not os.path.isfile(pdb2):
            ntf2 = tempfile.NamedTemporaryFile(mode='w')
            ntf2.write(pdb2)
            ntf2.flush()
            pdb2 = ntf2.name
        args = [pdb1, pdb2] + options.split()

        # mican calc
        process = subprocess.Popen([self.binary,'-z']+args,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        # output
        output, _ = process.communicate()
        outdict = output2dict(output.decode('utf-8'))

        # close tempfile
        if ntf1 is not None: ntf1.close()
        if ntf2 is not None: ntf2.close()
        # return
        return outdict if return_dict else Alignment(outdict)

