import pytest
import numpy as np
import pandas as pd

from mutationsPy.simulate_mut_matrix import sample_mut_mat_simple, sample_mut_mat_multinomial, gen_mut_matrix


def test_sample_mut_mat_simple():
    exposures = [[0.4, 0.6],
                [0.3, 0.7],
                [0.5, 0.5],
                [0.6, 0.4]]
    burdens = [5, 10, 20, 30]
    signatures = [[0.3, 0.3, 0.4], 
                  [0.1, 0.6, 0.3]]
    expected = np.array([[5*(0.3*0.4+0.1*0.6), 5*(0.3*0.4+0.6*0.6), 5*(0.4*0.4+0.3*0.6)],
                [10*(0.3*0.3+0.1*0.7), 10*(0.3*0.3+0.6*0.7), 10*(0.4*0.3+0.3*0.7)],
                [20*(0.3*0.5+0.1*0.5), 20*(0.3*0.5+0.6*0.5), 20*(0.4*0.5+0.3*0.5)],
                [30*(0.3*0.6+0.1*0.4), 30*(0.3*0.6+0.6*0.4), 30*(0.4*0.6+0.3*0.4)]])
    assert np.allclose(sample_mut_mat_simple(exposures, burdens, signatures), expected)
    
def test_sample_mut_mat_multinomial():
    exposures = [[0.4, 0.6],
                [0.3, 0.7],
                [0.5, 0.5],
                [0.6, 0.4]]
    burdens = [5, 10, 20, 30]
    signatures = [[0.3, 0.3, 0.4], 
                  [0.1, 0.6, 0.3]]
    
    np.random.seed(10)
    sample1_sig1 = np.random.multinomial(round(5*0.4), [0.3, 0.3, 0.4])
    sample1_sig2 = np.random.multinomial(round(5*0.6), [0.1, 0.6, 0.3])
    sample2_sig1 = np.random.multinomial(round(10*0.3), [0.3, 0.3, 0.4])
    sample2_sig2 = np.random.multinomial(round(10*0.7), [0.1, 0.6, 0.3])
    sample3_sig1 = np.random.multinomial(round(20*0.5), [0.3, 0.3, 0.4])
    sample3_sig2 = np.random.multinomial(round(20*0.5), [0.1, 0.6, 0.3])
    sample4_sig1 = np.random.multinomial(round(30*0.6), [0.3, 0.3, 0.4])
    sample4_sig2 = np.random.multinomial(round(30*0.4), [0.1, 0.6, 0.3])
    expected = np.array([
        sample1_sig1 + sample1_sig2, 
        sample2_sig1 + sample2_sig2, 
        sample3_sig1 + sample3_sig2, 
        sample4_sig1 + sample4_sig2
    ])
    np.random.seed(10)
    assert np.allclose(sample_mut_mat_multinomial(exposures, burdens, signatures), expected)
    
def test_gen_mut_matrix():
    exposure_matrix = pd.DataFrame(
        data = {'sample_id': [1, 2, 3, 4],
                'burdens': [5, 10, 20, 30],
                'sig1': [0.4, 0.3, 0.5, 0.6],
                'sig2': [0.6, 0.7, 0.5, 0.4]}
    )
    signatures = pd.DataFrame(
        data = {'Type': ['mut1', 'mut2', 'mut3'], 'sig1': [0.3, 0.3, 0.4], 'sig2': [0.1, 0.6, 0.3]}
    )
    assert list(gen_mut_matrix(exposure_matrix, signatures).columns) == ['sample_id', 'mut1', 'mut2', 'mut3']
    assert list(gen_mut_matrix(exposure_matrix, signatures, sample_func=sample_mut_mat_simple).columns) == ['sample_id', 'mut1', 'mut2', 'mut3']
    
def test_gen_mut_matrix_exposure_sum_unequal_one():
    exposure_matrix = pd.DataFrame(
        data = {'sample_id': [1, 2, 3, 4],
                'burdens': [5, 10, 20, 30],
                'sig1': [0.5, 0.3, 0.5, 0.6],
                'sig2': [0.6, 0.7, 0.5, 0.4]}
    )
    signatures = pd.DataFrame(
        data = {'Type': ['mut1', 'mut2', 'mut3'], 'sig1': [0.3, 0.3, 0.4], 'sig2': [0.1, 0.6, 0.3]}
    )
    with pytest.raises(ValueError):
        gen_mut_matrix(exposure_matrix, signatures)