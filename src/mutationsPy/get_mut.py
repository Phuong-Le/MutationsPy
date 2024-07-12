
def get_mut(context):
    """get the central mutation given the context 

    here mut refers to the base substitutions, eg C>G
    context refers to the substitution plus the flanking sequences, eg AT[C>G]GT
    seq refers to a sequence (a string), eg ATCGC
    Args:
        context (str): context of mutation

    Returns:
        str: the mutation
    """
    return context.split('[')[1].split(']')[0]

def get_context(seq, pos, ref, alt, kmer=3):
    """get the context given a reference sequence, the position and the state of the base at that position
    usually provided by a VCF format file 

    here mut refers to the base substitutions, eg C>G
    context refers to the substitution plus the flanking sequences, eg AT[C>G]GT
    seq refers to a sequence (a string), eg ATCGC
    Args:
        seq (str): reference sequence
        pos (int): position of the mutation
        ref (str): reference allele
        alt (str): mutated allele
        kmer (int, optional): kmer size. Defaults to 3.

    Raises:
        ValueError: kmer should be an odd number no less than 1
        KeyError: reference allele should match the base on the reference sequence

    Returns:
        str: sequence context of the mutation
    """
    if kmer % 2 != 1 or kmer < 1:
        raise ValueError("kmer should be an odd number no less than 1") 
    flank_size = kmer // 2
    if seq[pos] != ref:
        raise KeyError(f'position {pos} of the provided sequence is {seq[pos]},  not {ref}')
    if kmer == 1:
        return f'{ref}>{alt}'
    else:
        return f'{seq[(pos - flank_size):pos]}[{ref}>{alt}]{seq[(pos + 1):(pos + flank_size + 1)]}'

def get_wt_seq(mut):
    """take a mutation with sequence context and return the wildtype sequence

    Args:
        mut (str): mutation, for example A[C>G]T

    Returns:
        str: wildtype sequence, output for example above should be ACT
    """
    left, rest = mut.split("[") 
    central, right = rest.split("]")
    ref, alt = central.split(">")
    return left + ref + right 


def complementary_base(base):
    """get complementary base

    here mut refers to the base substitutions, eg C>G
    context refers to the substitution plus the flanking sequences, eg AT[C>G]GT
    seq refers to a sequence (a string), eg ATCGC
    Args:
        base (str): base of interest

    Returns:
        str: complementary base
    """
    complementary_bases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return complementary_bases[base]

def complementary_mut(mut):
    """get complementary mutation

    here mut refers to the base substitutions, eg C>G
    context refers to the substitution plus the flanking sequences, eg AT[C>G]GT
    seq refers to a sequence (a string), eg ATCGC
    Args:
        mut (str): mutation of interest

    Raises:
        ValueError: mutation has to be a substitution 

    Returns:
        str: complementary mutation 
    """
    if len(mut) != 3:
        raise ValueError("mut has to be a string of length 3")
    return f'{complementary_base(mut[0])}>{complementary_base(mut[2])}'
    
def complementary_seq(seq):
    """get complementary sequence

    here mut refers to the base substitutions, eg C>G
    context refers to the substitution plus the flanking sequences, eg AT[C>G]GT
    seq refers to a sequence (a string), eg ATCGC
    Args:
        seq (str): sequence of interest

    Returns:
        str: complementary sequence
    """
    rv_bases = [complementary_base(base) for base in seq[::-1]]
    return ''.join(rv_bases)

def rv_context(context):
    """get complementary context

    here mut refers to the base substitutions, eg C>G
    context refers to the substitution plus the flanking sequences, eg AT[C>G]GT
    seq refers to a sequence (a string), eg ATCGC
    Args:
        context (str): context of interest

    Returns:
        str: complementary mutation context
    """
    left, rest = context.split('[')
    mut, right = rest.split(']')
    return f'{complementary_seq(right)}[c{complementary_mut(mut)}]{complementary_seq(left)}'

def symmetrise_context(context):
    """get complementary context only if reference allele is either A or G

    here mut refers to the base substitutions, eg C>G
    context refers to the substitution plus the flanking sequences, eg AT[C>G]GT
    seq refers to a sequence (a string), eg ATCGC
    Args:
        context (str): context of interest

    Returns:
        str: "symmetrised" context
    """
    left, rest = context.split('[')
    mut, right = rest.split(']')
    result = context if mut[0] in ['C', 'T'] else f'{complementary_seq(right)}[{complementary_mut(mut)}]{complementary_seq(left)}'
    return result