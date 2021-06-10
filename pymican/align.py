import ctypes as ct
import os
import sys

dir_script = os.path.dirname(os.path.realpath(__file__))
LIBFILEPATH = os.path.abspath(dir_script+'/../mican.so')

libc = ct.cdll.LoadLibrary(LIBFILEPATH)
libc.main.restype = None
libc.main.argtypes = ct.c_int, ct.POINTER(ct.c_char_p)


class MICAN:
    def __init__(self):
        file1 = 'tests/test1.1.pdb'
        file2 = 'test/test1.2.pdb'
        stdout = os.dup(1)
        silent = os.open(os.devnull, os.O_WRONLY)
        os.dup2(silent, 1)
        args = (ct.c_char_p * 3)(b'mican',"{}".format(file1).encode(),b'tests/test1.2.pdb')
        libc.main(len(args), args)
        os.dup2(stdout, 1)


