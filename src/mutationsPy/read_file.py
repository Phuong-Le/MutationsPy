import pandas as pd

def read_vcf(vcf, pos_colname='POS', sep ='\t', **kwargs):
    """Read VCF format files, adjust the position so that it is 0-based indexing

    Args:
        vcf (path): path to VCF format file
        pos_colname (str, optional): name of the column that contains mutation positions. Defaults to 'POS'.
        sep (str, optional): delimiter of VCF file. Defaults to '\t'.

    Returns:
        a dataframe of VCF format with 0-based indexing
    """
    dt = pd.read_csv(vcf, sep = sep, **kwargs)
    dt[pos_colname] = dt[pos_colname] - 1
    return dt

def extract_seq_from_fasta(fasta, identifier):
    """get one sequence from fasta file based on the identifier

    Args:
        fasta (path): path to fasta file
        identifier (str): identifier of the sequence (without the ">")

    Returns:
        an uppercase sequence 
    """
    started = False 
    seqs = []
    identifier = '>' + identifier
    with open(fasta) as file:
        for line in file.readlines():
            if not line.strip(): # skip if line is empty
                continue
            if line.split()[0] == identifier:
                started = True
            if started == True:
                if line.startswith('>') and line.split()[0] != identifier:
                    break
                else:
                    if line.split()[0] != identifier:
                        seqs.append(line.strip('\n').strip())
    return ''.join(seqs).upper()