#script to repickle hmm_key.tsv for use in main.py after any alterations
import pickle
file = open('Databases_key.tsv', 'r')
dict1 = {}
for line in file:
    name = line.split('\t')[0]
    one = line.split('\t')[1]
    two = line.split('\t')[2]
    three = (line.split('\t')[3]).strip('\n')
    if not name in dict1:
        dict1[name] = {}
        dict1[name][one] = two +'\t' + three
    else:
        dict1[name][one] = two +'\t' + three
#print(dict1)
with open('Databases_key.pickle', 'wb') as handle:
    pickle.dump(dict1, handle, protocol=pickle.HIGHEST_PROTOCOL)
