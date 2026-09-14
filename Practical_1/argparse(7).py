import argparse
import sys
from Bio.PDB.PDBParser import PDBParser
import numpy as np

parser = argparse.ArgumentParser(
                                 prog='Exercise 7', 
                                 description='Print distances between all atom pairs of two given residues'
                                 )

# Use add_argument for all needed arguments
# https://docs.python.org/3/library/argparse.html

# Format for a optional text argument, default values can be indicated.
parser.add_argument('residue_1',
                    help='Name of the first residue'
                    )

parser.add_argument('residue_2',
                    help='Name of the second residue'
                    )

# Format for a required parameter, call will fail if empty. This only stores the file name, but specifiying type=file, the file is open and contents available. 
parser.add_argument('PDB_file',
                    help='Required PDB file for the program')

# Read command line into args

args = parser.parse_args()
       
# Print the parameters that has been read 
    
print ("\nSettings\n--------")

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

print ("\nSettings, again\n---------------")    

#print the variables once assigned
Res_1 = args.residue_1
Res_2 = args.residue_2
Pdb_file = args.PDB_file

print(Res_1, Res_2, Pdb_file) 

PDB_parser = PDBParser()

st = PDB_parser.get_structure(args.PDB_file, args.PDB_file)

def get_residue(structure, res_spec):
    """Parse 'CHAIN:NUM' or just 'NUM' and return the Biopython residue."""
    if ':' in res_spec:
        chain_id, resnum_str = res_spec.split(':')
    else:
        chain_id = next(iter(structure[0].child_dict))
        resnum_str = res_spec

    resnum = int(resnum_str)
    if chain_id not in structure[0]:
        sys.exit(f"Chain '{chain_id}' not found")

    chain = structure[0][chain_id]
    if resnum not in [res.id[1] for res in chain]:
        sys.exit(f"Residue {resnum} not found in chain {chain_id}")
    return chain[resnum]

res1 = get_residue(st, args.residue_1)
res2 = get_residue(st, args.residue_2)

print("Residue 1 is", res1.get_resname())
print("Residue 2 is", res2.get_resname())

print("\nAtom1 Atom2 dist1 dist2\n-------------------------")
pairs = []
for at1 in res1.get_atoms():      # Replace get_atoms with get_atom if you get an Error!
    for at2 in res2.get_atoms():
        dist = at2 - at1     # Direct procedure with (-) to compute distances
        vector = at2.coord - at1.coord  # Or using numpy coordinates
        distance = np.sqrt(np.sum(vector ** 2))
        pairs.append((at1, at2, dist, distance))

pairs.sort(key=lambda p: (p[0].get_serial_number(), p[1].get_serial_number()))

for at1, at2, dist, distance in pairs:
    print(at1, at2, dist, distance)

center = np.array([10, 10, 10])
print("\nDistance of res1 to {} \n".format(center))
res1_dists = []
for at1 in res1.get_atoms():
    vect = at1.coord - center
    distance = np.sqrt(np.sum(vect ** 2))
    res1_dists.append((at1, distance))

res1_dists.sort(key=lambda p: p[0].get_serial_number())

for at1, distance in res1_dists:
    print(at1, distance)

"""
Example usage:
1) Valid input: python 'argparse(7).py' 1 76 1ubq.pdb
Settings
--------
residue_1 : 1
residue_2 : 76
PDB_file  : 1ubq.pdb

Settings, again
---------------
1 76 1ubq.pdb
Residue 1 is MET
Residue 2 is GLY

Atom1 Atom2 dist1 dist2
-------------------------
<Atom N> <Atom N> 36.772102 36.772102
<Atom N> <Atom CA> 37.25674 37.25674
<Atom N> <Atom C> 38.474125 38.474125
<Atom N> <Atom O> 38.56514 38.56514
<Atom N> <Atom OXT> 39.289482 39.289482
<Atom CA> <Atom N> 36.630276 36.630276
<Atom CA> <Atom CA> 37.063484 37.063484
<Atom CA> <Atom C> 38.263958 38.263958
<Atom CA> <Atom O> 38.30952 38.30952
<Atom CA> <Atom OXT> 39.112232 39.112232
<Atom C> <Atom N> 35.3632 35.3632
<Atom C> <Atom CA> 35.772568 35.772568
<Atom C> <Atom C> 36.98727 36.987274
<Atom C> <Atom O> 37.0312 37.0312
<Atom C> <Atom OXT> 37.848557 37.848557
<Atom O> <Atom N> 34.419117 34.419117
<Atom O> <Atom CA> 34.858418 34.858418
<Atom O> <Atom C> 36.08381 36.08381
<Atom O> <Atom O> 36.155827 36.155827
<Atom O> <Atom OXT> 36.926052 36.926052
<Atom CB> <Atom N> 36.632633 36.632633
<Atom CB> <Atom CA> 37.063187 37.063187
<Atom CB> <Atom C> 38.224445 38.224445
<Atom CB> <Atom O> 38.238983 38.238983
<Atom CB> <Atom OXT> 39.07552 39.075516
<Atom CG> <Atom N> 35.322857 35.322857
<Atom CG> <Atom CA> 35.76603 35.76603
<Atom CG> <Atom C> 36.910564 36.910564
<Atom CG> <Atom O> 36.922817 36.922817
<Atom CG> <Atom OXT> 37.753517 37.753517
<Atom SD> <Atom N> 35.720596 35.720596
<Atom SD> <Atom CA> 36.165775 36.165775
<Atom SD> <Atom C> 37.258583 37.258583
<Atom SD> <Atom O> 37.2364 37.2364
<Atom SD> <Atom OXT> 38.098457 38.098457
<Atom CE> <Atom N> 34.135933 34.135933
<Atom CE> <Atom CA> 34.599823 34.599823
<Atom CE> <Atom C> 35.673862 35.673862
<Atom CE> <Atom O> 35.65397 35.65397
<Atom CE> <Atom OXT> 36.500217 36.500214

Distance of res1 to [10 10 10] 

<Atom N> 23.737175672088306
<Atom CA> 23.52403699425117
<Atom C> 24.591743612232424
<Atom O> 24.97703947233213
<Atom CB> 22.13870154645495
<Atom CG> 21.913744648642556
<Atom SD> 20.141395107391904
<Atom CE> 20.246788861024175

2) Valid input, residue + chain: 
python 'argparse(7).py' A:10 A:20 1ubq.pdb

Settings
--------
residue_1 : A:10
residue_2 : A:20
PDB_file  : 1ubq.pdb

Settings, again
---------------
A:10 A:20 1ubq.pdb
Residue 1 is GLY
Residue 2 is SER

Atom1 Atom2 dist1 dist2
-------------------------
<Atom N> <Atom N> 26.841194 26.841194
<Atom N> <Atom CA> 27.705828 27.705826
<Atom N> <Atom C> 27.01152 27.011518
<Atom N> <Atom O> 27.704405 27.704405
<Atom N> <Atom CB> 28.923367 28.923365
<Atom N> <Atom OG> 28.661741 28.661741
<Atom CA> <Atom N> 26.76021 26.76021
<Atom CA> <Atom CA> 27.68607 27.68607
<Atom CA> <Atom C> 27.0538 27.0538
<Atom CA> <Atom O> 27.799332 27.799332
<Atom CA> <Atom CB> 28.883638 28.883638
<Atom CA> <Atom OG> 28.57509 28.57509
<Atom C> <Atom N> 26.097872 26.097872
<Atom C> <Atom CA> 27.042477 27.042479
<Atom C> <Atom C> 26.412207 26.412207
<Atom C> <Atom O> 27.175257 27.175257
<Atom C> <Atom CB> 28.193848 28.19385
<Atom C> <Atom OG> 27.822264 27.822264
<Atom O> <Atom N> 26.35525 26.35525
<Atom O> <Atom CA> 27.34624 27.34624
<Atom O> <Atom C> 26.76351 26.76351
<Atom O> <Atom O> 27.56535 27.56535
<Atom O> <Atom CB> 28.469149 28.469149
<Atom O> <Atom OG> 28.049194 28.049194

Distance of res1 to [10 10 10] 

<Atom N> 39.46418688865108
<Atom CA> 39.183283863832415
<Atom C> 38.96855544105957
<Atom O> 39.107196877658936

3) Same residue: python 'argparse(7).py' 35 35 1ubq.pdb

Settings
--------
residue_1 : 35
residue_2 : 35
PDB_file  : 1ubq.pdb

Settings, again
---------------
35 35 1ubq.pdb
Residue 1 is GLY
Residue 2 is GLY

Atom1 Atom2 dist1 dist2
-------------------------
<Atom N> <Atom N> 0.0 0.0
<Atom N> <Atom CA> 1.4388051 1.4388051
<Atom N> <Atom C> 2.5115757 2.5115757
<Atom N> <Atom O> 3.6574357 3.6574357
<Atom CA> <Atom N> 1.4388051 1.4388051
<Atom CA> <Atom CA> 0.0 0.0
<Atom CA> <Atom C> 1.4781702 1.4781702
<Atom CA> <Atom O> 2.3734956 2.3734956
<Atom C> <Atom N> 2.5115757 2.5115757
<Atom C> <Atom CA> 1.4781702 1.4781702
<Atom C> <Atom C> 0.0 0.0
<Atom C> <Atom O> 1.2588807 1.2588807
<Atom O> <Atom N> 3.6574357 3.6574357
<Atom O> <Atom CA> 2.3734956 2.3734956
<Atom O> <Atom C> 1.2588807 1.2588807
<Atom O> <Atom O> 0.0 0.0

Distance of res1 to [10 10 10] 

<Atom N> 39.844092700591965
<Atom CA> 40.679445334042676
<Atom C> 40.42921512008421
<Atom O> 41.220835412978126

4) Invalid input (chain 999 doesn't exist): python 'argparse(7).py' 999 10 1ubq.pdb

Settings
--------
residue_1 : 999
residue_2 : 10
PDB_file  : 1ubq.pdb

Settings, again
---------------
999 10 1ubq.pdb
Residue 999 not found in chain A

"""