import argparse
import sys
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(
                                 prog='Exersice 1', 
                                 description='Determine the list of pairs of residues whose CA atoms are closer than a given distance'
                                 )

# Use add_argument for all needed arguments
# https://docs.python.org/3/library/argparse.html

# Format for a optional text argument, default values can be indicated.
parser.add_argument('distance', 
                    type=float,
                    help='Submit the distance between the CA atoms'
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
Variable1 = args.distance
Variable_Required = args.PDB_file

print(Variable1, Variable_Required) 

MAXDIST = Variable1  # Define distance for a  contact
if MAXDIST <= 0:
    sys.exit(f"Error: distance must be positive, got {MAXDIST}")

PDB_parser = PDBParser()

# load structure from PDB file

st = PDB_parser.get_structure(args.PDB_file, args.PDB_file)

select = []

#Select only CA atoms

for at in st.get_atoms():
    if at.id == 'CA':
        select.append(at)
        print(f"ATOM: {at.get_parent().get_resname()}, {at.get_parent().id[1]}, {at.id}")

# Preparing search
nbsearch = NeighborSearch(select)

print("NBSEARCH:")

contacts = [
    (a, b) if a.get_parent().id[1] <= b.get_parent().id[1] else (b, a)
    for a, b in nbsearch.search_all(MAXDIST)
]
contacts.sort(key=lambda pair: (pair[0].get_parent().id[1], pair[1].get_parent().id[1]))


ncontact = 1

for at1, at2 in contacts:
    print(f"Contact: {ncontact}")
    print(f"at1: {at1}, {at1.get_serial_number()}, {at1.get_parent().get_resname()}")
    print(f"at2: {at2}, {at2.get_serial_number()}, {at2.get_parent().get_resname()}")
    print()
    ncontact += 1

"""
Example usage: 
1) Boundary case: python 'argparse(1).py' 0 1ubq.pdb

Settings
--------
distance  : 0.0
PDB_file  : 1ubq.pdb

Settings, again
---------------
0.0 1ubq.pdb
Error: distance must be positive, got 0.0

2) Valid input: python 'argparse(1).py' 3.8 1ubq.pdb

Settings
--------
distance  : 3.8
PDB_file  : 1ubq.pdb

Settings, again
---------------
3.8 1ubq.pdb
ATOM: MET, 1, CA
ATOM: GLN, 2, CA
ATOM: ILE, 3, CA
ATOM: PHE, 4, CA
ATOM: VAL, 5, CA
ATOM: LYS, 6, CA
ATOM: THR, 7, CA
ATOM: LEU, 8, CA
ATOM: THR, 9, CA
ATOM: GLY, 10, CA
ATOM: LYS, 11, CA
ATOM: THR, 12, CA
ATOM: ILE, 13, CA
ATOM: THR, 14, CA
ATOM: LEU, 15, CA
ATOM: GLU, 16, CA
ATOM: VAL, 17, CA
ATOM: GLU, 18, CA
ATOM: PRO, 19, CA
ATOM: SER, 20, CA
ATOM: ASP, 21, CA
ATOM: THR, 22, CA
ATOM: ILE, 23, CA
ATOM: GLU, 24, CA
ATOM: ASN, 25, CA
ATOM: VAL, 26, CA
ATOM: LYS, 27, CA
ATOM: ALA, 28, CA
ATOM: LYS, 29, CA
ATOM: ILE, 30, CA
ATOM: GLN, 31, CA
ATOM: ASP, 32, CA
ATOM: LYS, 33, CA
ATOM: GLU, 34, CA
ATOM: GLY, 35, CA
ATOM: ILE, 36, CA
ATOM: PRO, 37, CA
ATOM: PRO, 38, CA
ATOM: ASP, 39, CA
ATOM: GLN, 40, CA
ATOM: GLN, 41, CA
ATOM: ARG, 42, CA
ATOM: LEU, 43, CA
ATOM: ILE, 44, CA
ATOM: PHE, 45, CA
ATOM: ALA, 46, CA
ATOM: GLY, 47, CA
ATOM: LYS, 48, CA
ATOM: GLN, 49, CA
ATOM: LEU, 50, CA
ATOM: GLU, 51, CA
ATOM: ASP, 52, CA
ATOM: GLY, 53, CA
ATOM: ARG, 54, CA
ATOM: THR, 55, CA
ATOM: LEU, 56, CA
ATOM: SER, 57, CA
ATOM: ASP, 58, CA
ATOM: TYR, 59, CA
ATOM: ASN, 60, CA
ATOM: ILE, 61, CA
ATOM: GLN, 62, CA
ATOM: LYS, 63, CA
ATOM: GLU, 64, CA
ATOM: SER, 65, CA
ATOM: THR, 66, CA
ATOM: LEU, 67, CA
ATOM: HIS, 68, CA
ATOM: LEU, 69, CA
ATOM: VAL, 70, CA
ATOM: LEU, 71, CA
ATOM: ARG, 72, CA
ATOM: LEU, 73, CA
ATOM: ARG, 74, CA
ATOM: GLY, 75, CA
ATOM: GLY, 76, CA
NBSEARCH:
Contact: 1
at1: <Atom CA>, 10, GLN
at2: <Atom CA>, 19, ILE

Contact: 2
at1: <Atom CA>, 69, THR
at2: <Atom CA>, 76, GLY

Contact: 3
at1: <Atom CA>, 76, GLY
at2: <Atom CA>, 80, LYS

Contact: 4
at1: <Atom CA>, 89, THR
at2: <Atom CA>, 96, ILE

Contact: 5
at1: <Atom CA>, 96, ILE
at2: <Atom CA>, 104, THR

Contact: 6
at1: <Atom CA>, 104, THR
at2: <Atom CA>, 111, LEU

Contact: 7
at1: <Atom CA>, 111, LEU
at2: <Atom CA>, 119, GLU

Contact: 8
at1: <Atom CA>, 119, GLU
at2: <Atom CA>, 128, VAL

Contact: 9
at1: <Atom CA>, 128, VAL
at2: <Atom CA>, 135, GLU

Contact: 10
at1: <Atom CA>, 144, PRO
at2: <Atom CA>, 151, SER

Contact: 11
at1: <Atom CA>, 151, SER
at2: <Atom CA>, 157, ASP

Contact: 12
at1: <Atom CA>, 157, ASP
at2: <Atom CA>, 165, THR

Contact: 13
at1: <Atom CA>, 165, THR
at2: <Atom CA>, 172, ILE

Contact: 14
at1: <Atom CA>, 189, ASN
at2: <Atom CA>, 197, VAL

Contact: 15
at1: <Atom CA>, 235, GLN
at2: <Atom CA>, 244, ASP

Contact: 16
at1: <Atom CA>, 261, GLU
at2: <Atom CA>, 270, GLY

Contact: 17
at1: <Atom CA>, 270, GLY
at2: <Atom CA>, 274, ILE

Contact: 18
at1: <Atom CA>, 349, PHE
at2: <Atom CA>, 360, ALA

Contact: 19
at1: <Atom CA>, 360, ALA
at2: <Atom CA>, 365, GLY

Contact: 20
at1: <Atom CA>, 448, ASP
at2: <Atom CA>, 456, TYR

Contact: 21
at1: <Atom CA>, 493, LYS
at2: <Atom CA>, 502, GLU

Contact: 22
at1: <Atom CA>, 502, GLU
at2: <Atom CA>, 511, SER

Contact: 23
at1: <Atom CA>, 557, LEU
at2: <Atom CA>, 565, ARG

Contact: 24
at1: <Atom CA>, 576, LEU
at2: <Atom CA>, 584, ARG

Contact: 25
at1: <Atom CA>, 584, ARG
at2: <Atom CA>, 595, GLY

Contact: 26
at1: <Atom CA>, 595, GLY
at2: <Atom CA>, 599, GLY






"""