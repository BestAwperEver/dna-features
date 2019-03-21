import pickle

with open("asdf.bin", "rb") as fp:
    prob_a = pickle.load(fp)

sorted_maxs_mins = [(332, 'min'), (568, 'min'), (1041, 'max'), (1202, 'min'), (1725, 'min')]

true_mins = []
true_maxs = []

ex_type = sorted_maxs_mins[0][1]
true_ex = sorted_maxs_mins[0][0]

for i in range(1, len(sorted_maxs_mins)):
    if sorted_maxs_mins[i][1] == ex_type:
        if ex_type == "max":
            if prob_a[sorted_maxs_mins[i][0]] > prob_a[true_ex]:
                true_ex = sorted_maxs_mins[i][0]
        else:
            if prob_a[sorted_maxs_mins[i][0]] < prob_a[true_ex]:
                true_ex = sorted_maxs_mins[i][0]
    else:
        if ex_type == "max":
            true_maxs.append(true_ex)
        else:
            true_mins.append(true_ex)
        ex_type = sorted_maxs_mins[i][1]
        true_ex = sorted_maxs_mins[i][0]