from itertools import combinations
from Levenshtein import distance as lev

def create_edgelist(cdr3, filename=None, method='HAMMING'):
    '''
    Create tab-separated edgelist of edges with HD = 1, from a set of sequences.    
    '''
    available = ['HAMMING', 'LEVENSHTEIN']
    assert method in available, f"Unknown distance metric, please choose one of the following:\n {available}"
    # Set makes sure there are no dupes
    cdr3 = set(cdr3)
    
    if method=='HAMMING':
        # Hashing
        cdr3hash = dict()
        for cdr in cdr3:
            for hash in (cdr[::2], cdr[1::2]):
                if hash not in cdr3hash:
                    cdr3hash[hash] = set()
                cdr3hash[hash].add(cdr)
                
        # Generate network
        edgelist = set()
        for hash in cdr3hash:
            if len(cdr3hash[hash]) >= 1:
                for cdr1 in cdr3hash[hash]:
                    for cdr2 in cdr3hash[hash]:
                        if cdr1 != cdr2:
                            if cdr1 <= cdr2:
                                if sum(ch1 != ch2 for ch1, ch2 in zip(cdr1, cdr2)) == 1:
                                    edgelist.add(cdr1 + "\t" + cdr2)
        
        # Save edgelist to file
        if filename is not None:
            with open(filename, 'w') as f:
                for edge in edgelist:
                    f.write('%s\n' % edge)
    
    elif method=='LEVENSHTEIN':
        edgelist = set()
        combos = [comb for comb in combinations(list(cdr3), 2)]
        for combo in combos:
            d = lev(combo[0],combo[1])
            if d == 1:
                edgelist.add(combo[0] + "\t" + combo[1])

    return edgelist