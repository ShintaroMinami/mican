import ctypes as ct

class PYALIGN(ct.Structure):
    _fields_ = [
        ('rot', (ct.c_float*4)*4),
        ('vec', ct.c_float*4),
        ('naa_align', ct.c_int),
        ('rmsd', ct.c_float),
        ('TMscore_ave', ct.c_float),
        ('TMscore_max', ct.c_float),
        ('TMscore_min', ct.c_float),
        ('seqID', ct.c_float)
        ]
