from itertools import product

def gen_context(kmer=3, symmetric=True):
    """Generate context given kmer size (ordered by the mutations then the flank sequence)

    Args:
        kmer (int, optional): kmer sizr. Defaults to 3.
        symmetric (bool, optional): whether mutation is symmetric. Defaults to True.

    Raises:
        ValueError: kmer has to be an odd number no less than 1

    Returns:
        a list of mutations 
    """
    if kmer % 2 != 1 or kmer < 1:
        raise ValueError("kmer should be an odd number no less than 1") 
    bases = 'ACGT'
    if symmetric:
        muts = [f'{ref}>{alt}' for ref, alt in product('CG', bases) if ref != alt]
    else:
        muts = [f'{ref}>{alt}' for ref, alt in product(bases, bases) if ref != alt]
    if kmer == 1: 
        return muts
    else:
        flank = [''.join(f) for f in product(bases, repeat=(kmer//2))]
        contexts = [f'{left}[{mut}]{right}' for mut, left, right in product(muts, flank, flank)]
        return contexts
        
        