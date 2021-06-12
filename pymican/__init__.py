import os
import subprocess
from typing import List, Union
from .parse_result import output2dict


dir_script = os.path.dirname(os.path.realpath(__file__))
BINFILEPATH = os.path.abspath(dir_script+'/bin/mican')

class mican:
    """
    MICAN: non-sequential alignment algorithm

    Attributes
    ----------
    binary : executable binary file path
    """
    def __init__(self, binary: str=BINFILEPATH):
        """
        Parameters
        ----------
        binary : executable binary file path
        """
        self.binary = binary
        return

    def align(self, pdb1: str, pdb2: str, options: Union[str, List[str]]=[])->dict:
        """
        Alignment calculation

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
        """
        # arguments
        options = options if type(options) == str else ' '.join(options)
        args = [pdb1, pdb2] + options.split()
        # mican calc
        process = subprocess.Popen([self.binary]+args,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        # output
        output, _ = process.communicate()
        output = output.decode('utf-8')
        # return
        return output2dict(output)
