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
	#	if power > 1:
	#		features[feature] = (power, p_ex - int_ex, p_in + int_in)
	#	else:
	#		features[feature] = (power, (p_ex + p_in)/2, (p_ex + p_in)/2)
	#else:
	#	power = (p_in - int_in) / (p_ex + int_ex)
	#	if power > 1:
	#		features[feature] = (power, p_ex + int_ex, p_in - int_in)
	#	else:
	#		features[feature] = (power, (p_ex + p_in)/2, (p_ex + p_in)/2)
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

f = open("files/transaction.txt")

e_i = float(f.readline())*200
i_e = float(f.readline())*200

f.close()

f = open("test_dna.txt")

dna = ""

for line in f:
	dna += line.strip().lower()

f.close()

dna_trans = [[1-e_i,e_i],[i_e,1-i_e]]

dna_states = ['E', 'I']

observations = ['A', 'C', 'G', 'T']
obs_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
states = ['H', 'L']


start_p = [.5, .5]
# start_p = {'H': 0.5, 'L': 0.5}
trans = [[.1, .9],[.1,.9]]
#trans = {
#    'H': {'H': 0.5, 'L': 0.5},
#    'L': {'H': 0.4, 'L': 0.6}
#}

emit = [[.2,.3,.3,.2],[.3,.2,.2,.3]]
#emit = {
#    'H': {'A': 0.2, 'C': 0.3, 'G': 0.3, 'T': 0.2},
#    'L': {'A': 0.3, 'C': 0.2, 'G': 0.2, 'T': 0.3}
#}

y = 'GGCACTGCGACTACGATACGAGCATATACTGACGACGACTGACGACGACTACTACTACGAAGCTAGCTACTAGCAGACATGACAGATGACTGACAGTCATGCATGACTGACTAA'
obs = [0 for i in range(len(y))]

for i in range(len(y)):
    obs[i] = obs_dict[y[i]]


def maxT_1A(T: np.array, i, A: np.array, j):
    k = T.shape[0]
    m = 0
    ind = 0
    for k0 in range(k):
        if T[k0][i] * A[k0][j] > m:
            m = T[k0][i] * A[k0][j]
            ind = k0
    return m, int(ind)


def argmax_z(T:np.array):
    k = T.shape[0]
    m = 0
    ind = 0
    for k1 in range(k):
        if T[k1][-1] > m:
            m = T[k1][-1]
            ind = k1
    return int(ind)


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

def viterbi_modified(obs: str, k: int, PI, A, B):
	# obs - наблюдения, k - количество скрытых состояний
	# PI - начальные вероятности
	# A - матрица перехода
	# B - матрица эмиссий
	T_1 = np.zeros((k, len(obs)))
	T_2 = np.zeros((k, len(obs)))
	z = [0 for i in range(len(obs))]
	# x = [0 for i in range(len(obs))]
	feature = most_powerful_feature(obs, N, 0)
	for i in range(k):
		# obs[0] заменить на наиболее мощный признак
		#T_1[i][0] = PI[i]*B[i][obs[0]]
		T_1[i][0] = PI[i]*B[feature][1+i]
		T_2[i][0] = 0
	for i in range(1, len(obs)):
		feature = most_powerful_feature(obs, N, i)
		for j in range(k):
			tup_max = maxT_1A(T_1, i-1, A, j)
			# obs[i] заменить на наиболее мощный признак
			T_1[j][i] = B[feature][1+j] * tup_max[0]
			T_2[j][i] = tup_max[1]
	z[-1] = argmax_z(T_1)
	# x_T = z[len(obs)] # вектор из 0 и 1
	for i in range(len(obs) - 1, 0, -1):
		z[i-1] = int(T_2[z[i]][i])
		# x[i-1] = [z[i-1]]  # вектор из 0 и 1
	return z

output = viterbi_modified(dna, 2, (.5, .5), dna_trans, features)

output = ''.join([dna_states[i] for i in output])

f = open("out.txt",'w')

f.write(output)

f.close()

def viterbi(obs, k, PI, A, B):
    T_1 = np.zeros((k, len(obs)))
    T_2 = np.zeros((k, len(obs)))
    z = [0 for i in range(len(obs))]
    # x = [0 for i in range(len(obs))]
    for i in range(k):
        T_1[i][0] = PI[i]*B[i][obs[0]]
        T_2[i][0] = 0
    for i in range(1, len(obs)):
        for j in range(k):
            tup_max = maxT_1A(T_1, i-1, A, j)
            T_1[j][i] = B[j][obs[i]] * tup_max[0]
            T_2[j][i] = tup_max[1]
    z[-1] = argmax_z(T_1)
    # x_T = z[len(obs)] # вектор из 0 и 1
    
    for i in range(len(obs) - 1, 0, -1):
        z[i-1] = int(T_2[z[i]][i])
        # x[i-1] = [z[i-1]]  # вектор из 0 и 1
    return z

#print(viterbi(obs, len(states), start_p, trans, emit))