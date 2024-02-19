import numpy as np
import pandas as pd 


def sample_mut_mat_simple(exposures, burdens, signatures):
    """calculate mutation matrix with exact original metrics

    Args:
        exposures (np.array convertible): exposures in proportion, where rows are the samples and columns are the signatures
        burdens (np.array convertible): a vector of burdens for each sample
        signatures (np.array convertible): signatures with rows are the signatures, and columns are the proportion of mutation classes

    Returns:
        _np.array_: mutation matrix with dimension compatible with SigProfilerExtractor
    """
    burden_by_signatures = np.array(exposures).T * np.array(burdens) 
    mut_mat = np.matmul(burden_by_signatures.T, signatures)
    return mut_mat

def sample_mut_mat_multinomial(exposures, burdens, signatures):
    """calculate mutation matrix with exact original metrics

    Args:
        exposures (np.array convertible): exposures in proportion, where rows are the samples and columns are the signatures
        burdens (np.array convertible): a vector of burdens for each sample
        signatures (np.array convertible): signatures with rows are the signatures, and columns are the proportion of mutation classes

    Returns:
        _np.array_: mutation matrix with dimension compatible with SigProfilerExtractor
    """
    exposures = np.array(exposures)
    mut_mat = []
    for i, burden in enumerate(burdens): 
        exposure = exposures[i]
        burden_by_sig = [np.random.multinomial(round(burden*exposure[j]), signatures[j]) for j, _ in enumerate(signatures)]
        burden_by_sig = np.array(burden_by_sig).sum(axis=0)
        mut_mat.append(burden_by_sig)
    return np.array(mut_mat)

def gen_mut_matrix(exposure_matrix, signatures, sample_func = sample_mut_mat_multinomial):
    """_summary_

    Args:
        exposure_matrix (pd.DataFrame): columns include 'burdens' (mutation burdens), 'sample_id' and 'signatures' (in proportion)
        signatures (pd.DataFrame): downloaded from COSMIC
        sample_func (function, optional): _description_. Defaults to sample_mut_mat_multinomial.

    Returns:
        _type_: _description_
    """
    burdens = exposure_matrix['burdens']
    mut_class = list(signatures['Type'])
    present_signatures = list(set(exposure_matrix.columns) & set(signatures.columns))
    exposures = exposure_matrix[present_signatures]
    if not np.allclose(exposures.sum(axis = 1), 1):
        raise ValueError('exposures from at least one sample do not add up to 1')
    signatures = np.array(signatures[present_signatures]).T
    mut_mat = sample_func(exposures, burdens, signatures)
    mut_mat = pd.DataFrame(mut_mat, columns = mut_class)
    mut_mat['sample_id'] = exposure_matrix['sample_id']
    mut_mat = mut_mat[['sample_id'] + mut_class]
    return round(mut_mat)