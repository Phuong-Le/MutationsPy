
import pytest
import pandas as pd
from tempfile import NamedTemporaryFile

from mutationsPy.mut_matrix import sigprofiler_to_hdp_no_rearrangement, sigprofiler_to_hdp, hdp_to_sigprofiler, concat_sigprofiler_mutmats

import os
import sys
current_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(current_dir)

def test_sigprofiler_to_hdp_no_rearrangement():
    sigprofiler_path = 'tests/test_data/mut_matrix/sigprofiler_mutmat.txt'
    with NamedTemporaryFile(delete=False, mode = 'w+t') as hdp_outpath:    
        sigprofiler_to_hdp_no_rearrangement(sigprofiler_path=sigprofiler_path, hdp_outpath=hdp_outpath.name)
        result = pd.read_csv(hdp_outpath.name, sep = '\t')
        print(result)
        expected = pd.read_csv('tests/test_data/mut_matrix/hdp_mutmat_no_rearrangement.txt', sep = '\t')
        print(expected)
        assert result.equals(expected)


def test_sigprofiler_to_hdp():
    sigprofiler_path = 'tests/test_data/mut_matrix/sigprofiler_mutmat.txt'
    with NamedTemporaryFile(delete=False, mode = 'w+t') as hdp_outpath:    
        sigprofiler_to_hdp(sigprofiler_path=sigprofiler_path, hdp_outpath=hdp_outpath.name)
        result = pd.read_csv(hdp_outpath.name, sep = '\t')
        print(result)
        expected = pd.read_csv('tests/test_data/mut_matrix/hdp_mutmat.txt', sep = '\t')
        print(expected)
        assert result.equals(expected)
    
def test_hdp_to_sigprofiler():
    hdp_path = 'tests/test_data/mut_matrix/hdp_mutmat.txt'
    with NamedTemporaryFile(delete=False, mode = 'w+t') as outpath:    
        hdp_to_sigprofiler(hdp_path=hdp_path, sigprofiler_outpath=outpath.name)
        result = pd.read_csv(outpath.name, sep = '\t')
        print(result)
        expected = pd.read_csv('tests/test_data/mut_matrix/sigprofiler_mutmat.txt', sep = '\t')
        print(expected)
        assert result.equals(expected)
