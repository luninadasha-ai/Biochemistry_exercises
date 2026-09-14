import argparse
import sys
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(
                                 prog='Exercise 5', 
                                 description='DGenerate a list of backbone connectivity (i.e. which residues are linked by ordinary peptide bonds).'
                                 )

# Use add_argument for all needed arguments
# https://docs.python.org/3/library/argparse.html

# Format for a optional text argument, default values can be indicated.
parser.add_argument(
    '-d', '--distance',
    type=float,
    default=2.5,
    help='Distance for a peptide, default - 2.5'
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
Dist = args.distance
Pdb_file = args.PDB_file

print(Dist, Pdb_file) 

MAXDIST = Dist  # Define distance for a  contact

PDB_parser = PDBParser()

# load structure from PDB file

st = PDB_parser.get_structure(args.PDB_file, args.PDB_file)

select = []

elems = ['C', 'N']

for at in st.get_atoms():
    if at.id in elems:
        select.append(at)
        print(f"ATOM: {at.get_parent().get_resname()}, {at.get_parent().id[1]}, {at.id}")

# Preparing search
nbsearch = NeighborSearch(select)

print("NBSEARCH:")

contacts = []

for a, b in nbsearch.search_all(MAXDIST):
    if a.id == 'C' and b.id == 'N':
        c_atom, n_atom = a, b
    elif a.id == 'N' and b.id == 'C':
        c_atom, n_atom = b, a
    else:
        continue  # not a C-N pair

    res_c = c_atom.get_parent()
    res_n = n_atom.get_parent()

    if res_c == res_n:
        continue  # same residue, not a peptide bond

    if res_c.get_parent().id != res_n.get_parent().id:
        continue  # different chains, not a peptide bond

    contacts.append((c_atom, n_atom))

contacts.sort(key=lambda pair: pair[0].get_parent().id[1])

ncontact = 1

for c_atom, n_atom in contacts:
    res_c = c_atom.get_parent()
    res_n = n_atom.get_parent()
    print(f"Bond {ncontact}: "
          f"{res_c.get_resname()} {res_c.id[1]} (C)  <-->  "
          f"{res_n.get_resname()} {res_n.id[1]} (N)")
    ncontact += 1

"""
Example usage: 

1) Normal run: python 'argparse(5).py' 1ubq.pdb
Settings
--------
distance  : 2.5
PDB_file  : 1ubq.pdb

Settings, again
---------------
2.5 1ubq.pdb
ATOM: MET, 1, N
ATOM: MET, 1, C
ATOM: GLN, 2, N
ATOM: GLN, 2, C
ATOM: ILE, 3, N
ATOM: ILE, 3, C
ATOM: PHE, 4, N
ATOM: PHE, 4, C
ATOM: VAL, 5, N
ATOM: VAL, 5, C
ATOM: LYS, 6, N
ATOM: LYS, 6, C
ATOM: THR, 7, N
ATOM: THR, 7, C
ATOM: LEU, 8, N
ATOM: LEU, 8, C
ATOM: THR, 9, N
ATOM: THR, 9, C
ATOM: GLY, 10, N
ATOM: GLY, 10, C
ATOM: LYS, 11, N
ATOM: LYS, 11, C
ATOM: THR, 12, N
ATOM: THR, 12, C
ATOM: ILE, 13, N
ATOM: ILE, 13, C
ATOM: THR, 14, N
ATOM: THR, 14, C
ATOM: LEU, 15, N
ATOM: LEU, 15, C
ATOM: GLU, 16, N
ATOM: GLU, 16, C
ATOM: VAL, 17, N
ATOM: VAL, 17, C
ATOM: GLU, 18, N
ATOM: GLU, 18, C
ATOM: PRO, 19, N
ATOM: PRO, 19, C
ATOM: SER, 20, N
ATOM: SER, 20, C
ATOM: ASP, 21, N
ATOM: ASP, 21, C
ATOM: THR, 22, N
ATOM: THR, 22, C
ATOM: ILE, 23, N
ATOM: ILE, 23, C
ATOM: GLU, 24, N
ATOM: GLU, 24, C
ATOM: ASN, 25, N
ATOM: ASN, 25, C
ATOM: VAL, 26, N
ATOM: VAL, 26, C
ATOM: LYS, 27, N
ATOM: LYS, 27, C
ATOM: ALA, 28, N
ATOM: ALA, 28, C
ATOM: LYS, 29, N
ATOM: LYS, 29, C
ATOM: ILE, 30, N
ATOM: ILE, 30, C
ATOM: GLN, 31, N
ATOM: GLN, 31, C
ATOM: ASP, 32, N
ATOM: ASP, 32, C
ATOM: LYS, 33, N
ATOM: LYS, 33, C
ATOM: GLU, 34, N
ATOM: GLU, 34, C
ATOM: GLY, 35, N
ATOM: GLY, 35, C
ATOM: ILE, 36, N
ATOM: ILE, 36, C
ATOM: PRO, 37, N
ATOM: PRO, 37, C
ATOM: PRO, 38, N
ATOM: PRO, 38, C
ATOM: ASP, 39, N
ATOM: ASP, 39, C
ATOM: GLN, 40, N
ATOM: GLN, 40, C
ATOM: GLN, 41, N
ATOM: GLN, 41, C
ATOM: ARG, 42, N
ATOM: ARG, 42, C
ATOM: LEU, 43, N
ATOM: LEU, 43, C
ATOM: ILE, 44, N
ATOM: ILE, 44, C
ATOM: PHE, 45, N
ATOM: PHE, 45, C
ATOM: ALA, 46, N
ATOM: ALA, 46, C
ATOM: GLY, 47, N
ATOM: GLY, 47, C
ATOM: LYS, 48, N
ATOM: LYS, 48, C
ATOM: GLN, 49, N
ATOM: GLN, 49, C
ATOM: LEU, 50, N
ATOM: LEU, 50, C
ATOM: GLU, 51, N
ATOM: GLU, 51, C
ATOM: ASP, 52, N
ATOM: ASP, 52, C
ATOM: GLY, 53, N
ATOM: GLY, 53, C
ATOM: ARG, 54, N
ATOM: ARG, 54, C
ATOM: THR, 55, N
ATOM: THR, 55, C
ATOM: LEU, 56, N
ATOM: LEU, 56, C
ATOM: SER, 57, N
ATOM: SER, 57, C
ATOM: ASP, 58, N
ATOM: ASP, 58, C
ATOM: TYR, 59, N
ATOM: TYR, 59, C
ATOM: ASN, 60, N
ATOM: ASN, 60, C
ATOM: ILE, 61, N
ATOM: ILE, 61, C
ATOM: GLN, 62, N
ATOM: GLN, 62, C
ATOM: LYS, 63, N
ATOM: LYS, 63, C
ATOM: GLU, 64, N
ATOM: GLU, 64, C
ATOM: SER, 65, N
ATOM: SER, 65, C
ATOM: THR, 66, N
ATOM: THR, 66, C
ATOM: LEU, 67, N
ATOM: LEU, 67, C
ATOM: HIS, 68, N
ATOM: HIS, 68, C
ATOM: LEU, 69, N
ATOM: LEU, 69, C
ATOM: VAL, 70, N
ATOM: VAL, 70, C
ATOM: LEU, 71, N
ATOM: LEU, 71, C
ATOM: ARG, 72, N
ATOM: ARG, 72, C
ATOM: LEU, 73, N
ATOM: LEU, 73, C
ATOM: ARG, 74, N
ATOM: ARG, 74, C
ATOM: GLY, 75, N
ATOM: GLY, 75, C
ATOM: GLY, 76, N
ATOM: GLY, 76, C
NBSEARCH:
Bond 1: MET 1 (C)  <-->  GLN 2 (N)
Bond 2: GLN 2 (C)  <-->  ILE 3 (N)
Bond 3: ILE 3 (C)  <-->  PHE 4 (N)
Bond 4: PHE 4 (C)  <-->  VAL 5 (N)
Bond 5: VAL 5 (C)  <-->  LYS 6 (N)
Bond 6: LYS 6 (C)  <-->  THR 7 (N)
Bond 7: THR 7 (C)  <-->  LEU 8 (N)
Bond 8: LEU 8 (C)  <-->  THR 9 (N)
Bond 9: THR 9 (C)  <-->  GLY 10 (N)
Bond 10: GLY 10 (C)  <-->  LYS 11 (N)
Bond 11: LYS 11 (C)  <-->  THR 12 (N)
Bond 12: THR 12 (C)  <-->  ILE 13 (N)
Bond 13: ILE 13 (C)  <-->  THR 14 (N)
Bond 14: THR 14 (C)  <-->  LEU 15 (N)
Bond 15: LEU 15 (C)  <-->  GLU 16 (N)
Bond 16: GLU 16 (C)  <-->  VAL 17 (N)
Bond 17: VAL 17 (C)  <-->  GLU 18 (N)
Bond 18: GLU 18 (C)  <-->  PRO 19 (N)
Bond 19: PRO 19 (C)  <-->  SER 20 (N)
Bond 20: SER 20 (C)  <-->  ASP 21 (N)
Bond 21: ASP 21 (C)  <-->  THR 22 (N)
Bond 22: THR 22 (C)  <-->  ILE 23 (N)
Bond 23: ILE 23 (C)  <-->  GLU 24 (N)
Bond 24: GLU 24 (C)  <-->  ASN 25 (N)
Bond 25: ASN 25 (C)  <-->  VAL 26 (N)
Bond 26: VAL 26 (C)  <-->  LYS 27 (N)
Bond 27: LYS 27 (C)  <-->  ALA 28 (N)
Bond 28: ALA 28 (C)  <-->  LYS 29 (N)
Bond 29: LYS 29 (C)  <-->  ILE 30 (N)
Bond 30: ILE 30 (C)  <-->  GLN 31 (N)
Bond 31: GLN 31 (C)  <-->  ASP 32 (N)
Bond 32: ASP 32 (C)  <-->  LYS 33 (N)
Bond 33: LYS 33 (C)  <-->  GLU 34 (N)
Bond 34: GLU 34 (C)  <-->  GLY 35 (N)
Bond 35: GLY 35 (C)  <-->  ILE 36 (N)
Bond 36: ILE 36 (C)  <-->  PRO 37 (N)
Bond 37: PRO 37 (C)  <-->  PRO 38 (N)
Bond 38: PRO 38 (C)  <-->  ASP 39 (N)
Bond 39: ASP 39 (C)  <-->  GLN 40 (N)
Bond 40: GLN 40 (C)  <-->  GLN 41 (N)
Bond 41: GLN 41 (C)  <-->  ARG 42 (N)
Bond 42: ARG 42 (C)  <-->  LEU 43 (N)
Bond 43: LEU 43 (C)  <-->  ILE 44 (N)
Bond 44: ILE 44 (C)  <-->  PHE 45 (N)
Bond 45: PHE 45 (C)  <-->  ALA 46 (N)
Bond 46: ALA 46 (C)  <-->  GLY 47 (N)
Bond 47: GLY 47 (C)  <-->  LYS 48 (N)
Bond 48: LYS 48 (C)  <-->  GLN 49 (N)
Bond 49: GLN 49 (C)  <-->  LEU 50 (N)
Bond 50: LEU 50 (C)  <-->  GLU 51 (N)
Bond 51: GLU 51 (C)  <-->  ASP 52 (N)
Bond 52: ASP 52 (C)  <-->  GLY 53 (N)
Bond 53: GLY 53 (C)  <-->  ARG 54 (N)
Bond 54: ARG 54 (C)  <-->  THR 55 (N)
Bond 55: THR 55 (C)  <-->  LEU 56 (N)
Bond 56: LEU 56 (C)  <-->  SER 57 (N)
Bond 57: SER 57 (C)  <-->  ASP 58 (N)
Bond 58: ASP 58 (C)  <-->  TYR 59 (N)
Bond 59: TYR 59 (C)  <-->  ASN 60 (N)
Bond 60: ASN 60 (C)  <-->  ILE 61 (N)
Bond 61: ILE 61 (C)  <-->  GLN 62 (N)
Bond 62: GLN 62 (C)  <-->  LYS 63 (N)
Bond 63: LYS 63 (C)  <-->  GLU 64 (N)
Bond 64: GLU 64 (C)  <-->  SER 65 (N)
Bond 65: SER 65 (C)  <-->  THR 66 (N)
Bond 66: THR 66 (C)  <-->  LEU 67 (N)
Bond 67: LEU 67 (C)  <-->  HIS 68 (N)
Bond 68: HIS 68 (C)  <-->  LEU 69 (N)
Bond 69: LEU 69 (C)  <-->  VAL 70 (N)
Bond 70: VAL 70 (C)  <-->  LEU 71 (N)
Bond 71: LEU 71 (C)  <-->  ARG 72 (N)
Bond 72: ARG 72 (C)  <-->  LEU 73 (N)
Bond 73: LEU 73 (C)  <-->  ARG 74 (N)
Bond 74: ARG 74 (C)  <-->  GLY 75 (N)
Bond 75: GLY 75 (C)  <-->  GLY 76 (N)

2) With a different cut off distance: 'argparse(5).py' 1ubq.pdb -d 2.0
Settings
--------
distance  : 2.0
PDB_file  : 1ubq.pdb

Settings, again
---------------
2.0 1ubq.pdb
ATOM: MET, 1, N
ATOM: MET, 1, C
ATOM: GLN, 2, N
ATOM: GLN, 2, C
ATOM: ILE, 3, N
ATOM: ILE, 3, C
ATOM: PHE, 4, N
ATOM: PHE, 4, C
ATOM: VAL, 5, N
ATOM: VAL, 5, C
ATOM: LYS, 6, N
ATOM: LYS, 6, C
ATOM: THR, 7, N
ATOM: THR, 7, C
ATOM: LEU, 8, N
ATOM: LEU, 8, C
ATOM: THR, 9, N
ATOM: THR, 9, C
ATOM: GLY, 10, N
ATOM: GLY, 10, C
ATOM: LYS, 11, N
ATOM: LYS, 11, C
ATOM: THR, 12, N
ATOM: THR, 12, C
ATOM: ILE, 13, N
ATOM: ILE, 13, C
ATOM: THR, 14, N
ATOM: THR, 14, C
ATOM: LEU, 15, N
ATOM: LEU, 15, C
ATOM: GLU, 16, N
ATOM: GLU, 16, C
ATOM: VAL, 17, N
ATOM: VAL, 17, C
ATOM: GLU, 18, N
ATOM: GLU, 18, C
ATOM: PRO, 19, N
ATOM: PRO, 19, C
ATOM: SER, 20, N
ATOM: SER, 20, C
ATOM: ASP, 21, N
ATOM: ASP, 21, C
ATOM: THR, 22, N
ATOM: THR, 22, C
ATOM: ILE, 23, N
ATOM: ILE, 23, C
ATOM: GLU, 24, N
ATOM: GLU, 24, C
ATOM: ASN, 25, N
ATOM: ASN, 25, C
ATOM: VAL, 26, N
ATOM: VAL, 26, C
ATOM: LYS, 27, N
ATOM: LYS, 27, C
ATOM: ALA, 28, N
ATOM: ALA, 28, C
ATOM: LYS, 29, N
ATOM: LYS, 29, C
ATOM: ILE, 30, N
ATOM: ILE, 30, C
ATOM: GLN, 31, N
ATOM: GLN, 31, C
ATOM: ASP, 32, N
ATOM: ASP, 32, C
ATOM: LYS, 33, N
ATOM: LYS, 33, C
ATOM: GLU, 34, N
ATOM: GLU, 34, C
ATOM: GLY, 35, N
ATOM: GLY, 35, C
ATOM: ILE, 36, N
ATOM: ILE, 36, C
ATOM: PRO, 37, N
ATOM: PRO, 37, C
ATOM: PRO, 38, N
ATOM: PRO, 38, C
ATOM: ASP, 39, N
ATOM: ASP, 39, C
ATOM: GLN, 40, N
ATOM: GLN, 40, C
ATOM: GLN, 41, N
ATOM: GLN, 41, C
ATOM: ARG, 42, N
ATOM: ARG, 42, C
ATOM: LEU, 43, N
ATOM: LEU, 43, C
ATOM: ILE, 44, N
ATOM: ILE, 44, C
ATOM: PHE, 45, N
ATOM: PHE, 45, C
ATOM: ALA, 46, N
ATOM: ALA, 46, C
ATOM: GLY, 47, N
ATOM: GLY, 47, C
ATOM: LYS, 48, N
ATOM: LYS, 48, C
ATOM: GLN, 49, N
ATOM: GLN, 49, C
ATOM: LEU, 50, N
ATOM: LEU, 50, C
ATOM: GLU, 51, N
ATOM: GLU, 51, C
ATOM: ASP, 52, N
ATOM: ASP, 52, C
ATOM: GLY, 53, N
ATOM: GLY, 53, C
ATOM: ARG, 54, N
ATOM: ARG, 54, C
ATOM: THR, 55, N
ATOM: THR, 55, C
ATOM: LEU, 56, N
ATOM: LEU, 56, C
ATOM: SER, 57, N
ATOM: SER, 57, C
ATOM: ASP, 58, N
ATOM: ASP, 58, C
ATOM: TYR, 59, N
ATOM: TYR, 59, C
ATOM: ASN, 60, N
ATOM: ASN, 60, C
ATOM: ILE, 61, N
ATOM: ILE, 61, C
ATOM: GLN, 62, N
ATOM: GLN, 62, C
ATOM: LYS, 63, N
ATOM: LYS, 63, C
ATOM: GLU, 64, N
ATOM: GLU, 64, C
ATOM: SER, 65, N
ATOM: SER, 65, C
ATOM: THR, 66, N
ATOM: THR, 66, C
ATOM: LEU, 67, N
ATOM: LEU, 67, C
ATOM: HIS, 68, N
ATOM: HIS, 68, C
ATOM: LEU, 69, N
ATOM: LEU, 69, C
ATOM: VAL, 70, N
ATOM: VAL, 70, C
ATOM: LEU, 71, N
ATOM: LEU, 71, C
ATOM: ARG, 72, N
ATOM: ARG, 72, C
ATOM: LEU, 73, N
ATOM: LEU, 73, C
ATOM: ARG, 74, N
ATOM: ARG, 74, C
ATOM: GLY, 75, N
ATOM: GLY, 75, C
ATOM: GLY, 76, N
ATOM: GLY, 76, C
NBSEARCH:
Bond 1: MET 1 (C)  <-->  GLN 2 (N)
Bond 2: GLN 2 (C)  <-->  ILE 3 (N)
Bond 3: ILE 3 (C)  <-->  PHE 4 (N)
Bond 4: PHE 4 (C)  <-->  VAL 5 (N)
Bond 5: VAL 5 (C)  <-->  LYS 6 (N)
Bond 6: LYS 6 (C)  <-->  THR 7 (N)
Bond 7: THR 7 (C)  <-->  LEU 8 (N)
Bond 8: LEU 8 (C)  <-->  THR 9 (N)
Bond 9: THR 9 (C)  <-->  GLY 10 (N)
Bond 10: GLY 10 (C)  <-->  LYS 11 (N)
Bond 11: LYS 11 (C)  <-->  THR 12 (N)
Bond 12: THR 12 (C)  <-->  ILE 13 (N)
Bond 13: ILE 13 (C)  <-->  THR 14 (N)
Bond 14: THR 14 (C)  <-->  LEU 15 (N)
Bond 15: LEU 15 (C)  <-->  GLU 16 (N)
Bond 16: GLU 16 (C)  <-->  VAL 17 (N)
Bond 17: VAL 17 (C)  <-->  GLU 18 (N)
Bond 18: GLU 18 (C)  <-->  PRO 19 (N)
Bond 19: PRO 19 (C)  <-->  SER 20 (N)
Bond 20: SER 20 (C)  <-->  ASP 21 (N)
Bond 21: ASP 21 (C)  <-->  THR 22 (N)
Bond 22: THR 22 (C)  <-->  ILE 23 (N)
Bond 23: ILE 23 (C)  <-->  GLU 24 (N)
Bond 24: GLU 24 (C)  <-->  ASN 25 (N)
Bond 25: ASN 25 (C)  <-->  VAL 26 (N)
Bond 26: VAL 26 (C)  <-->  LYS 27 (N)
Bond 27: LYS 27 (C)  <-->  ALA 28 (N)
Bond 28: ALA 28 (C)  <-->  LYS 29 (N)
Bond 29: LYS 29 (C)  <-->  ILE 30 (N)
Bond 30: ILE 30 (C)  <-->  GLN 31 (N)
Bond 31: GLN 31 (C)  <-->  ASP 32 (N)
Bond 32: ASP 32 (C)  <-->  LYS 33 (N)
Bond 33: LYS 33 (C)  <-->  GLU 34 (N)
Bond 34: GLU 34 (C)  <-->  GLY 35 (N)
Bond 35: GLY 35 (C)  <-->  ILE 36 (N)
Bond 36: ILE 36 (C)  <-->  PRO 37 (N)
Bond 37: PRO 37 (C)  <-->  PRO 38 (N)
Bond 38: PRO 38 (C)  <-->  ASP 39 (N)
Bond 39: ASP 39 (C)  <-->  GLN 40 (N)
Bond 40: GLN 40 (C)  <-->  GLN 41 (N)
Bond 41: GLN 41 (C)  <-->  ARG 42 (N)
Bond 42: ARG 42 (C)  <-->  LEU 43 (N)
Bond 43: LEU 43 (C)  <-->  ILE 44 (N)
Bond 44: ILE 44 (C)  <-->  PHE 45 (N)
Bond 45: PHE 45 (C)  <-->  ALA 46 (N)
Bond 46: ALA 46 (C)  <-->  GLY 47 (N)
Bond 47: GLY 47 (C)  <-->  LYS 48 (N)
Bond 48: LYS 48 (C)  <-->  GLN 49 (N)
Bond 49: GLN 49 (C)  <-->  LEU 50 (N)
Bond 50: LEU 50 (C)  <-->  GLU 51 (N)
Bond 51: GLU 51 (C)  <-->  ASP 52 (N)
Bond 52: ASP 52 (C)  <-->  GLY 53 (N)
Bond 53: GLY 53 (C)  <-->  ARG 54 (N)
Bond 54: ARG 54 (C)  <-->  THR 55 (N)
Bond 55: THR 55 (C)  <-->  LEU 56 (N)
Bond 56: LEU 56 (C)  <-->  SER 57 (N)
Bond 57: SER 57 (C)  <-->  ASP 58 (N)
Bond 58: ASP 58 (C)  <-->  TYR 59 (N)
Bond 59: TYR 59 (C)  <-->  ASN 60 (N)
Bond 60: ASN 60 (C)  <-->  ILE 61 (N)
Bond 61: ILE 61 (C)  <-->  GLN 62 (N)
Bond 62: GLN 62 (C)  <-->  LYS 63 (N)
Bond 63: LYS 63 (C)  <-->  GLU 64 (N)
Bond 64: GLU 64 (C)  <-->  SER 65 (N)
Bond 65: SER 65 (C)  <-->  THR 66 (N)
Bond 66: THR 66 (C)  <-->  LEU 67 (N)
Bond 67: LEU 67 (C)  <-->  HIS 68 (N)
Bond 68: HIS 68 (C)  <-->  LEU 69 (N)
Bond 69: LEU 69 (C)  <-->  VAL 70 (N)
Bond 70: VAL 70 (C)  <-->  LEU 71 (N)
Bond 71: LEU 71 (C)  <-->  ARG 72 (N)
Bond 72: ARG 72 (C)  <-->  LEU 73 (N)
Bond 73: LEU 73 (C)  <-->  ARG 74 (N)
Bond 74: ARG 74 (C)  <-->  GLY 75 (N)
Bond 75: GLY 75 (C)  <-->  GLY 76 (N)

3) Wrong input: python 'argparse(5).py' 1ubq.pdb -d 'not_a_number'
usage: Exercise 5 [-h] [-d DISTANCE] PDB_file
Exercise 5: error: argument -d/--distance: invalid float value: 'not_a_number'





"""
