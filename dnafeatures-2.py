import numpy as np
from collections import OrderedDict as OD
from scipy import stats

N = 7
prob = .93

f = open("files/features.deep-"+str(N)+".txt")

z = stats.norm.ppf(prob)

features = {}

for line in f:
	raw_feature = line.strip('\n').split()
	feature = raw_feature[0]
	p_ex = float(raw_feature[1])
	int_ex = z*float(raw_feature[2])
	p_in = float(raw_feature[3])
	int_in = z*float(raw_feature[4])
	power = 0
	if p_ex > p_in:
		power = (p_ex - int_ex) / (p_in + int_in)
		if power > 1:
			features[feature] = (power, power / (1 + power), 1 / (1 + power))
		else:
			features[feature] = (power, .5, .5)
	else:
		power = (p_in - int_in) / (p_ex + int_ex)
		if power > 1:
			features[feature] = (power, 1 / (1 + power), power / (1 + power))
		else:
			features[feature] = (power, .5, .5)

f.close()

f = open("files/important.features.deep-7.prob-90.0.threshold-2.5.txt")

important_features = {}

for line in f:
	raw_feature = line.strip('\n').split()
	feature = raw_feature[0]
	p_ex = float(raw_feature[1])
	int_ex = z*float(raw_feature[2])
	p_in = float(raw_feature[3])
	int_in = z*float(raw_feature[4])
	power = 0
	if p_ex > p_in:
		power = (p_ex - int_ex) / (p_in + int_in)
		important_features[feature] = (power, 'e')
	else:
		power = (p_in - int_in) / (p_ex + int_ex)
		important_features[feature] = (power, 'i')

f = open("files/transaction.txt")

e_i = float(f.readline())
i_e = float(f.readline())

f.close()

dna_trans = [[1-e_i,e_i],[i_e,1-i_e]]

def most_powerful_feature(obs: str, N: int, index: int) -> str:
	i = 1 # max(5, len(obs) - index)
	pow = 0
	feature = obs[index]
	while i <= N and index + i <= len(obs):
		if features[obs[index:index+i]][0] > pow:
			pow = features[obs[index:index+i]][0] 
			feature = obs[index:index+i]
		i += 1
	return feature

f = open("test_dna.txt")

type = "exon"
dna = ""
for line in f:
	l = line.strip()
	if len(l) > 1:
		changed = False
		if l[0].islower():
			if type == "exon":
				changed = True
			type = "intron"
		else:
			if type == "intron":
				changed = True
			type = "exon"
		if changed:
			imp_features = {'e': 0, 'i': 0}
			for i in range(len(dna)):
				feat = most_powerful_feature(dna, N, i)
				if important_features.get(feat) != None:
					imp_features[important_features[feat][1]] += 1
			if type == "exon":
				type1 = "intron"
			else:
				type1 = "exon"
			print(imp_features, type1, (imp_features['i'] < imp_features['e']) == (type1 == "exon"), sep='\t')
			dna = ""
		dna += l.lower()


f.close()