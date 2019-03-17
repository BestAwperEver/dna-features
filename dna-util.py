import numpy as np
import itertools

def str2subset(s, l):
    return [s[i:i+l] for i in range(len(s)-l+1)]


def str2subsets(s, min_l, max_l):
    assert max_l <= len(s)
    sets = {}
    for i in range(min_l, max_l+1):
        sets[i] = str2subset(s, i)
    return list(itertools.chain.from_iterable(list(sets.values())))


def get_doc(f, upper):
    s = ''
    c = f.read(1)
    while c:
        if c not in 'acgtACGT':
            c = f.read(1)
            continue
        if c.isupper() != upper:
            return s, c
        s += c
        c = f.read(1)
    return s, c


def make_docs(filePath):
    exons = []
    introns = []
    f = open(filePath)
    c = f.read(1)
    while (c):
        upper = c.isupper()
        pair = get_doc(f, upper)
        doc = c + pair[0]
        if upper:
            exons.append(doc)
        else:
            introns.append(doc)
        c = pair[1]
    return exons, introns


def get_subset_docs(filePath, min_l, max_l):
    exons, introns = make_docs(filePath)
    exons = [str2subsets(s, min_l, max_l) for s in exons]
    introns = [str2subsets(s, min_l, max_l) for s in introns]
    return exons, introns


docs = make_docs("files/dna.txt")
print([str2subset(s, 2) for s in docs[1]])
