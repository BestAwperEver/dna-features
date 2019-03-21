import csv
from nltk.stem.porter import PorterStemmer
from math import log
import gensim.parsing.preprocessing as prep
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.metrics import f1_score
import numpy as np
import itertools
from sklearn.naive_bayes import MultinomialNB
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.pipeline import Pipeline


def str2subset(s, l):
    return [s[i:i+l] for i in range(len(s)-l+1)]


def str2subsets_dict(s, min_l, max_l):
    assert max_l <= len(s)
    sets = {}
    for i in range(min_l, max_l+1):
        sets[i] = str2subset(s, i)
#     return list(itertools.chain.from_iterable(list(sets.values())))
    return sets


def str2subsets_str(s, min_l, max_l):
    assert max_l <= len(s)
    sets = []
    for i in range(min_l, max_l+1):
        sets += str2subset(s, i)
    return ' '.join(sets)


def get_docs(f, upper, min_len = 50):
    s = ''
    c = f.read(1)
    while c:
        if c not in 'acgtACGT':
            c = f.read(1)
            continue
        if c.isupper() != upper:
            break
        s += c
        c = f.read(1)
    docs = []
    k = 0
    while k+min_len < len(s):
        docs.append(s[k:k+min_len])
        k += min_len
    if len(docs) > 0:
        docs[-1] += s[k:-1]
    else:
        docs.append(s)
    return docs, c


def make_docs(filePath):
    exons = []
    introns = []
    f = open(filePath)
    c = f.read(1)
    while (c):
        upper = c.isupper()
        pair = get_docs(f, upper, 100)
        doc = c + pair[0][0]
        if upper:
            exons.append(doc)
        else:
            introns.append(doc)
        for i in range(1, len(pair[0])):
            doc = pair[0][i]
            if upper:
                exons.append(doc)
            else:
                introns.append(doc)
        c = pair[1]
    return exons, introns


def get_subset_docs(filePath, min_l, max_l):
    exons, introns = make_docs(filePath)
    exons = [str2subsets_str(s, min_l, max_l) for s in exons]
    introns = [str2subsets_str(s, min_l, max_l) for s in introns]
    return exons, introns


def get_entropy(p):
    assert p >=0 and p <= 1
    if p < 10e-10 or p > 1-10e-10:
        return 0
    return -p*log(p)-(1-p)*log(1-p)


def get_sum_entropy(clf, dnas):
    sum_len = sum([len(s) for s in dnas])
    probs = clf.predict_proba(dnas).T[0].tolist()
    h = 0
    for i in range(len(dnas)):
        h += get_entropy(probs[i])*float(len(dnas[i]))/sum_len
    return h, probs


def optimal_split(s:str, clf):
    ig = 0
    h0 = get_entropy(clf.predict_proba([s])[0][0])
    i0 = 0
    probs0 = []
    for i in range(len(s)):
        h, probs = get_sum_entropy(clf, [s[:i], s[i:]])
        h /= len(s)
        if h0 - h > ig:
            ig = h0 - h
            i0 = i
            probs0 = probs
    return s[:i0], s[i0:], probs0


def split_rec(s:str, clf):
    s1, s2, probs = optimal_split(s, clf)
    if len(s1) == 0:
        return [s2]
    return split_rec(s1) + split_rec(s2)


exons, introns = get_subset_docs(r"C:\Users\admin\PycharmProjects\dna\files\dna.txt", 1, 7)
# data = pd.DataFrame(data=np.array([exons+introns,
#                           np.concatenate((np.ones_like(exons, dtype=np.int8),
#                                           np.zeros_like(introns, dtype=np.int8))).tolist()]).T,
#                     columns=['dna', 'exon'])
#
# text_clf = Pipeline([('vect', CountVectorizer()),
#                      ('tfidf', TfidfTransformer()),
#                      ('clf', MultinomialNB()),
# ])
#
# kf = KFold(n_splits=10, shuffle=True)
# f1_scores = []
# for train, test in kf.split(data):
#     X = data['dna'][train]
#     Y = data['exon'][train]
#     text_clf.fit(data['dna'], data['exon'])
#     predicted = text_clf.predict(data['dna'][test])
#     f1_scores.append(f1_score(data['exon'][test], predicted, average='micro'))
# print("F-score: ", sum(f1_scores)/len(f1_scores))

important_features = []
with open(r"C:\Users\admin\PycharmProjects\dna\files\important.features.deep-7.m-100.prob-find-0.999.prob-72.0.threshold-3.0.txt", 'r') as f:
    for line in f:
        if line[0] != 'e':
            important_features.append(line.split()[0])

exons_only_important = []
for s in exons:
    s_imp = []
    for ss in s.split():
        if ss in important_features:
            s_imp.append(ss)
    exons_only_important.append(' '.join(s_imp))

introns_only_important = []
for s in introns:
    s_imp = []
    for ss in s:
        if ss in important_features:
            s_imp.append(ss)
    introns_only_important.append(s_imp)

test_dna = "agataaaaatagtaaagatattcatatttatacagctttacaagttgaaacatcctttca\
tttatgaagaattaaaaggggtaccctttttagagaaaaggagagcatgtaaacttcgag\
gaaattgatatgtataattttataaaacagggcttgcgcttttttttttttgagacagag\
tttcgctcttgttgcccaggctggagtgcaatggtgcaacctcggctcaccgcaacctcc\
tcctcccgagttcaagtgattctcctgcctcagcctgctgaatagctgggattacaggca\
tgtgccaccacacctggctacttttgtgttttttttacttttatatattttttttttgtt\
tagtagagacagggtttctccattttggtcaggctggtcttgaactcccgacctcagatg\
atctgcccgcctcagcctcccaaagtgctgggattacaggcgtgagccactgtgcctggc\
caggggttgtgctttttaaatttcaattttatttttgctaagtatttattctttgatag\
TTCATGTTTTGGGAAAAGAACAGGCTTCACCTAAAAACGTAAAAATGGAAATTGGTAAAA\
CTGAAACTTTTTCTGATGTTCCTGTGAAAACAAATATAGAAGTTTGTTCTACTTACTCCA\
AAGATTCAGAAAACTACTTTGAAACAGAAGCAGTAGAAATTGCTAAAGCTTTTATGGAAG\
ATGATGAACTGACAGATTCTAAACTGCCAAGTCATGCCACACATTCTCTTTTTACATGTC\
CCGAAAATGAGGAAATGGTTTTGTCAAATTCAAGAATTGGAAAAAGAAGAGGAGAGCCCC\
TTATCTTAGTGGGAGAACCCTCAATCAAAAGAAACTTATTAAATGAATTTGACAGGATAATAGAAAATCAAG\
AAAAATCCTTAAAGGCTTCAAAAAGCACTCCAGATG\
GCACAATAAAAGATCGAAGATTGTTTATGCATCATGTTTCTTTAGAGCCGATTACCTGTG\
TACCCTTTCG"


text_clf = Pipeline([('vect', CountVectorizer()),
                     ('tfidf', TfidfTransformer()),
                     ('clf', MultinomialNB()),
])

data_imp = pd.DataFrame(data=np.array([exons_only_important+introns_only_important,
                          np.concatenate((np.ones_like(exons_only_important, dtype=np.int8),
                                          np.zeros_like(introns_only_important, dtype=np.int8))).tolist()]).T,
                    columns=['dna', 'exon'])