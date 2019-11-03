import json

from molecules import collect_molecules

mol = collect_molecules(cg=False)

for key in mol:
    mol_dict = {}
    mol_dict.update({'head' : mol[key][0]})
    tails = []
    for lst in mol[key][1]:
        tails.append(lst.tolist())
    mol_dict.update({'tails' : tails})
    mol_dict.update({'n_atoms' : mol[key][2]})
    mol_dict.update({'name' : key})
    mol[key] = mol_dict

for key in mol:
    print(key)
    for attr in mol[key]:
        print('\t' + str(attr) + '\t' + str(mol[key][attr]))

for key in mol:
    with open('molecules/atomistic/{}.json'.format(key), 'w') as f:
        json.dump(mol[key], f,  indent=4, sort_keys=True)
