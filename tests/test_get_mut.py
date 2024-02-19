import pytest
from mutationsPy.get_mut import get_mut, get_context, complementary_seq, symmetrise_context

def test_get_mut():
    assert get_mut('A[T>C]G') == 'T>C'
    assert get_mut('GAT[G>A]GCT') == 'G>A'
    
def test_get_context():
    assert get_context('ACGTAC', pos=2, ref='G', alt='T') == 'C[G>T]T'
    with pytest.raises(ValueError):
        get_context('ACGTAC', pos=2, ref='G', alt='T', kmer=0)
    with pytest.raises(ValueError):
        get_context('ACGTAC', pos=2, ref='G', alt='T', kmer=4)
    with pytest.raises(KeyError):
        get_context('ACGTAC', pos=2, ref='T', alt='T')
        
def test_complementary_seq():
    assert complementary_seq('ACGGT') == 'ACCGT'
    
def test_symmetrise_context():
    assert symmetrise_context('AGG[A>C]CGT') == 'ACG[T>G]CCT'
    assert symmetrise_context('ACG[T>G]CCT') == 'ACG[T>G]CCT'