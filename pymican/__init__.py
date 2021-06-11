import ctypes as ct
import os
from typing import List
import numpy as np
from .cstruct import PYALIGN

dir_script = os.path.dirname(os.path.realpath(__file__))
LIBFILEPATH = os.path.abspath(dir_script+'/../mican.so')

class mican:
    def __init__(self):
        # mican setup
        self.libc = ct.cdll.LoadLibrary(LIBFILEPATH)
        self.libc.pymain.restype = PYALIGN
        self.libc.pymain.argtypes = ct.c_int, ct.POINTER(ct.c_char_p)
        return

    def align(self, file1: str, file2: str, options: List[str]=[])->dict:
        # arguments
        args = [a.encode() for a in ['mican', file1, file2] + options]
        args = (ct.c_char_p * len(args))(*args)
        # mican calc
        output = self.libc.pymain(len(args), args)
        # prepare output dict
        outdict = {
            'naa_align': output.naa_align,
            'rmsd': output.rmsd,
            'TMscore_ave': output.TMscore_ave,
            'TMscore_max': output.TMscore_max,
            'TMscore_min': output.TMscore_min,
            'seqID': output.seqID,
            'rot': np.ctypeslib.as_array(output.rot, shape=(4,4))[1:,1:],
            'vec': np.ctypeslib.as_array(output.vec, shape=(4))[1:]
        }
        # output
        return outdict
