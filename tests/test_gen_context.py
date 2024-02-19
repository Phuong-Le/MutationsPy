from mutationsPy.gen_context import gen_context

def test_gen_context():
    assert gen_context(1) == ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    assert len(gen_context(3)) == 6 * (4**(3 - 1))
    assert len(gen_context(3, symmetric=False)) == 12 * (4**(3 - 1))
    assert len(gen_context(5)) == 6 * (4**(5 - 1))
    assert len(gen_context(5, symmetric=False)) == 12 * (4**(5 - 1))