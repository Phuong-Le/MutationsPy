import pytest
from tempfile import NamedTemporaryFile
import pandas as pd

from mutationsPy.read_file import read_vcf, extract_seq_from_fasta

def test_read_vcf():
    with NamedTemporaryFile(delete=False, mode = 'w+t') as vcf:
        vcf.writelines("CHROM\tPOS\tREF\tALT\nchr1\t30\tA\tC\nchr15\t10\tG\tC\n")
        vcf.seek(0)
        expected = pd.DataFrame({'CHROM': ['chr1', 'chr15'], 'POS': [29, 9], 'REF': ['A', 'G'], 'ALT': ['C', 'C']})
        assert read_vcf(vcf.name).equals(expected)
    
def test_read_vcf_with_kwargs():
    with NamedTemporaryFile(delete=False, mode = 'w+t') as vcf:
        vcf.writelines("CHROM\tPOS\tREF\tALT\tFILTER\nchr1\t30\tA\tC\tPASS\nchr15\t10\tG\tC\tNot Passed\n")
        vcf.seek(0)
        expected = pd.DataFrame({'CHROM': ['chr1', 'chr15'], 'POS': [29, 9], 'REF': ['A', 'G'], 'ALT': ['C', 'C']})
        assert read_vcf(vcf.name, usecols = ['CHROM', 'POS', 'REF', 'ALT']).equals(expected)

def test_extract_seq_from_fasta():
    with NamedTemporaryFile(delete=False, mode = 'w+t') as fasta:
        fasta.writelines(">chr1 \nACGTCA  \nAGCTACTAGCATANN\n>chr10\nGCTACNNNacgtaA\n\n>chr2   \nGCTACTGCA")
        fasta.seek(0)
        # test for line concatenation
        assert extract_seq_from_fasta(fasta.name, identifier='chr1') == 'ACGTCAAGCTACTAGCATANN'
        # test for uppercase and skip empty lines
        assert extract_seq_from_fasta(fasta.name, identifier='chr10') == 'GCTACNNNACGTAA'
        assert extract_seq_from_fasta(fasta.name, identifier='chr2') == 'GCTACTGCA'
        