import argparse
import sys
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(
                                 prog='Exercise 3', 
                                 description='Determine all possible hydrogen bonds (Polar atoms at less than 3.5 Å)'
                                 )

# Use add_argument for all needed arguments
# https://docs.python.org/3/library/argparse.html

# Format for a optional text argument, default values can be indicated.
parser.add_argument(
    '-d', '--distance',
    type=float,
    default=3.5,
    help='Distance for a hydrogen bond, default - 3.5'
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

polar_elems = ['O', 'N', 'S']

for at in st.get_atoms():
    if at.element in polar_elems: 
        select.append(at)
        print(f"ATOM: {at.get_parent().get_resname()}, {at.get_parent().id[1]}, {at.id}")

# Preparing search
nbsearch = NeighborSearch(select)

print("NBSEARCH:")

contacts = [
    (a, b) if a.get_parent().id[1] <= b.get_parent().id[1] else (b, a)
    for a, b in nbsearch.search_all(MAXDIST)
    if a.get_parent() != b.get_parent()
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

1) Default valid input: python 'argparse(3).py' 1ubq.pdb
Matches all 328 contacts, so I will not paste it here.

2) Different distance: python 'argparse(3).py' 1ubq.pdb -d 2.5

Settings
--------
distance  : 2.5
PDB_file  : 1ubq.pdb

Settings, again
---------------
2.5 1ubq.pdb
ATOM: MET, 1, N
ATOM: MET, 1, O
ATOM: MET, 1, SD
ATOM: GLN, 2, N
ATOM: GLN, 2, O
ATOM: GLN, 2, OE1
ATOM: GLN, 2, NE2
ATOM: ILE, 3, N
ATOM: ILE, 3, O
ATOM: PHE, 4, N
ATOM: PHE, 4, O
ATOM: VAL, 5, N
ATOM: VAL, 5, O
ATOM: LYS, 6, N
ATOM: LYS, 6, O
ATOM: LYS, 6, NZ
ATOM: THR, 7, N
ATOM: THR, 7, O
ATOM: THR, 7, OG1
ATOM: LEU, 8, N
ATOM: LEU, 8, O
ATOM: THR, 9, N
ATOM: THR, 9, O
ATOM: THR, 9, OG1
ATOM: GLY, 10, N
ATOM: GLY, 10, O
ATOM: LYS, 11, N
ATOM: LYS, 11, O
ATOM: LYS, 11, NZ
ATOM: THR, 12, N
ATOM: THR, 12, O
ATOM: THR, 12, OG1
ATOM: ILE, 13, N
ATOM: ILE, 13, O
ATOM: THR, 14, N
ATOM: THR, 14, O
ATOM: THR, 14, OG1
ATOM: LEU, 15, N
ATOM: LEU, 15, O
ATOM: GLU, 16, N
ATOM: GLU, 16, O
ATOM: GLU, 16, OE1
ATOM: GLU, 16, OE2
ATOM: VAL, 17, N
ATOM: VAL, 17, O
ATOM: GLU, 18, N
ATOM: GLU, 18, O
ATOM: GLU, 18, OE1
ATOM: GLU, 18, OE2
ATOM: PRO, 19, N
ATOM: PRO, 19, O
ATOM: SER, 20, N
ATOM: SER, 20, O
ATOM: SER, 20, OG
ATOM: ASP, 21, N
ATOM: ASP, 21, O
ATOM: ASP, 21, OD1
ATOM: ASP, 21, OD2
ATOM: THR, 22, N
ATOM: THR, 22, O
ATOM: THR, 22, OG1
ATOM: ILE, 23, N
ATOM: ILE, 23, O
ATOM: GLU, 24, N
ATOM: GLU, 24, O
ATOM: GLU, 24, OE1
ATOM: GLU, 24, OE2
ATOM: ASN, 25, N
ATOM: ASN, 25, O
ATOM: ASN, 25, OD1
ATOM: ASN, 25, ND2
ATOM: VAL, 26, N
ATOM: VAL, 26, O
ATOM: LYS, 27, N
ATOM: LYS, 27, O
ATOM: LYS, 27, NZ
ATOM: ALA, 28, N
ATOM: ALA, 28, O
ATOM: LYS, 29, N
ATOM: LYS, 29, O
ATOM: LYS, 29, NZ
ATOM: ILE, 30, N
ATOM: ILE, 30, O
ATOM: GLN, 31, N
ATOM: GLN, 31, O
ATOM: GLN, 31, OE1
ATOM: GLN, 31, NE2
ATOM: ASP, 32, N
ATOM: ASP, 32, O
ATOM: ASP, 32, OD1
ATOM: ASP, 32, OD2
ATOM: LYS, 33, N
ATOM: LYS, 33, O
ATOM: LYS, 33, NZ
ATOM: GLU, 34, N
ATOM: GLU, 34, O
ATOM: GLU, 34, OE1
ATOM: GLU, 34, OE2
ATOM: GLY, 35, N
ATOM: GLY, 35, O
ATOM: ILE, 36, N
ATOM: ILE, 36, O
ATOM: PRO, 37, N
ATOM: PRO, 37, O
ATOM: PRO, 38, N
ATOM: PRO, 38, O
ATOM: ASP, 39, N
ATOM: ASP, 39, O
ATOM: ASP, 39, OD1
ATOM: ASP, 39, OD2
ATOM: GLN, 40, N
ATOM: GLN, 40, O
ATOM: GLN, 40, OE1
ATOM: GLN, 40, NE2
ATOM: GLN, 41, N
ATOM: GLN, 41, O
ATOM: GLN, 41, OE1
ATOM: GLN, 41, NE2
ATOM: ARG, 42, N
ATOM: ARG, 42, O
ATOM: ARG, 42, NE
ATOM: ARG, 42, NH1
ATOM: ARG, 42, NH2
ATOM: LEU, 43, N
ATOM: LEU, 43, O
ATOM: ILE, 44, N
ATOM: ILE, 44, O
ATOM: PHE, 45, N
ATOM: PHE, 45, O
ATOM: ALA, 46, N
ATOM: ALA, 46, O
ATOM: GLY, 47, N
ATOM: GLY, 47, O
ATOM: LYS, 48, N
ATOM: LYS, 48, O
ATOM: LYS, 48, NZ
ATOM: GLN, 49, N
ATOM: GLN, 49, O
ATOM: GLN, 49, OE1
ATOM: GLN, 49, NE2
ATOM: LEU, 50, N
ATOM: LEU, 50, O
ATOM: GLU, 51, N
ATOM: GLU, 51, O
ATOM: GLU, 51, OE1
ATOM: GLU, 51, OE2
ATOM: ASP, 52, N
ATOM: ASP, 52, O
ATOM: ASP, 52, OD1
ATOM: ASP, 52, OD2
ATOM: GLY, 53, N
ATOM: GLY, 53, O
ATOM: ARG, 54, N
ATOM: ARG, 54, O
ATOM: ARG, 54, NE
ATOM: ARG, 54, NH1
ATOM: ARG, 54, NH2
ATOM: THR, 55, N
ATOM: THR, 55, O
ATOM: THR, 55, OG1
ATOM: LEU, 56, N
ATOM: LEU, 56, O
ATOM: SER, 57, N
ATOM: SER, 57, O
ATOM: SER, 57, OG
ATOM: ASP, 58, N
ATOM: ASP, 58, O
ATOM: ASP, 58, OD1
ATOM: ASP, 58, OD2
ATOM: TYR, 59, N
ATOM: TYR, 59, O
ATOM: TYR, 59, OH
ATOM: ASN, 60, N
ATOM: ASN, 60, O
ATOM: ASN, 60, OD1
ATOM: ASN, 60, ND2
ATOM: ILE, 61, N
ATOM: ILE, 61, O
ATOM: GLN, 62, N
ATOM: GLN, 62, O
ATOM: GLN, 62, OE1
ATOM: GLN, 62, NE2
ATOM: LYS, 63, N
ATOM: LYS, 63, O
ATOM: LYS, 63, NZ
ATOM: GLU, 64, N
ATOM: GLU, 64, O
ATOM: GLU, 64, OE1
ATOM: GLU, 64, OE2
ATOM: SER, 65, N
ATOM: SER, 65, O
ATOM: SER, 65, OG
ATOM: THR, 66, N
ATOM: THR, 66, O
ATOM: THR, 66, OG1
ATOM: LEU, 67, N
ATOM: LEU, 67, O
ATOM: HIS, 68, N
ATOM: HIS, 68, O
ATOM: HIS, 68, ND1
ATOM: HIS, 68, NE2
ATOM: LEU, 69, N
ATOM: LEU, 69, O
ATOM: VAL, 70, N
ATOM: VAL, 70, O
ATOM: LEU, 71, N
ATOM: LEU, 71, O
ATOM: ARG, 72, N
ATOM: ARG, 72, O
ATOM: ARG, 72, NE
ATOM: ARG, 72, NH1
ATOM: ARG, 72, NH2
ATOM: LEU, 73, N
ATOM: LEU, 73, O
ATOM: ARG, 74, N
ATOM: ARG, 74, O
ATOM: ARG, 74, NE
ATOM: ARG, 74, NH1
ATOM: ARG, 74, NH2
ATOM: GLY, 75, N
ATOM: GLY, 75, O
ATOM: GLY, 76, N
ATOM: GLY, 76, O
ATOM: GLY, 76, OXT
ATOM: HOH, 77, O
ATOM: HOH, 78, O
ATOM: HOH, 79, O
ATOM: HOH, 80, O
ATOM: HOH, 81, O
ATOM: HOH, 82, O
ATOM: HOH, 83, O
ATOM: HOH, 84, O
ATOM: HOH, 85, O
ATOM: HOH, 86, O
ATOM: HOH, 87, O
ATOM: HOH, 88, O
ATOM: HOH, 89, O
ATOM: HOH, 90, O
ATOM: HOH, 91, O
ATOM: HOH, 92, O
ATOM: HOH, 93, O
ATOM: HOH, 94, O
ATOM: HOH, 95, O
ATOM: HOH, 96, O
ATOM: HOH, 97, O
ATOM: HOH, 98, O
ATOM: HOH, 99, O
ATOM: HOH, 100, O
ATOM: HOH, 101, O
ATOM: HOH, 102, O
ATOM: HOH, 103, O
ATOM: HOH, 104, O
ATOM: HOH, 105, O
ATOM: HOH, 106, O
ATOM: HOH, 107, O
ATOM: HOH, 108, O
ATOM: HOH, 109, O
ATOM: HOH, 110, O
ATOM: HOH, 111, O
ATOM: HOH, 112, O
ATOM: HOH, 113, O
ATOM: HOH, 114, O
ATOM: HOH, 115, O
ATOM: HOH, 116, O
ATOM: HOH, 117, O
ATOM: HOH, 118, O
ATOM: HOH, 119, O
ATOM: HOH, 120, O
ATOM: HOH, 121, O
ATOM: HOH, 122, O
ATOM: HOH, 123, O
ATOM: HOH, 124, O
ATOM: HOH, 125, O
ATOM: HOH, 126, O
ATOM: HOH, 127, O
ATOM: HOH, 128, O
ATOM: HOH, 129, O
ATOM: HOH, 130, O
ATOM: HOH, 131, O
ATOM: HOH, 132, O
ATOM: HOH, 133, O
ATOM: HOH, 134, O
NBSEARCH:
Contact: 1
at1: <Atom O>, 4, MET
at2: <Atom N>, 9, GLN

Contact: 2
at1: <Atom N>, 1, MET
at2: <Atom O>, 637, HOH

Contact: 3
at1: <Atom O>, 12, GLN
at2: <Atom N>, 18, ILE

Contact: 4
at1: <Atom O>, 21, ILE
at2: <Atom N>, 26, PHE

Contact: 5
at1: <Atom O>, 29, PHE
at2: <Atom N>, 37, VAL

Contact: 6
at1: <Atom O>, 40, VAL
at2: <Atom N>, 44, LYS

Contact: 7
at1: <Atom O>, 47, LYS
at2: <Atom N>, 53, THR

Contact: 8
at1: <Atom O>, 56, THR
at2: <Atom N>, 60, LEU

Contact: 9
at1: <Atom O>, 63, LEU
at2: <Atom N>, 68, THR

Contact: 10
at1: <Atom O>, 71, THR
at2: <Atom N>, 75, GLY

Contact: 11
at1: <Atom O>, 78, GLY
at2: <Atom N>, 79, LYS

Contact: 12
at1: <Atom O>, 82, LYS
at2: <Atom N>, 88, THR

Contact: 13
at1: <Atom O>, 91, THR
at2: <Atom N>, 95, ILE

Contact: 14
at1: <Atom O>, 98, ILE
at2: <Atom N>, 103, THR

Contact: 15
at1: <Atom O>, 106, THR
at2: <Atom N>, 110, LEU

Contact: 16
at1: <Atom O>, 113, LEU
at2: <Atom N>, 118, GLU

Contact: 17
at1: <Atom O>, 121, GLU
at2: <Atom N>, 127, VAL

Contact: 18
at1: <Atom O>, 130, VAL
at2: <Atom N>, 134, GLU

Contact: 19
at1: <Atom O>, 137, GLU
at2: <Atom N>, 143, PRO

Contact: 20
at1: <Atom OE1>, 141, GLU
at2: <Atom O>, 625, HOH

Contact: 21
at1: <Atom O>, 146, PRO
at2: <Atom N>, 150, SER

Contact: 22
at1: <Atom O>, 153, SER
at2: <Atom N>, 156, ASP

Contact: 23
at1: <Atom O>, 159, ASP
at2: <Atom N>, 164, THR

Contact: 24
at1: <Atom O>, 167, THR
at2: <Atom N>, 171, ILE

Contact: 25
at1: <Atom O>, 174, ILE
at2: <Atom N>, 179, GLU

Contact: 26
at1: <Atom O>, 182, GLU
at2: <Atom N>, 188, ASN

Contact: 27
at1: <Atom O>, 191, ASN
at2: <Atom N>, 196, VAL

Contact: 28
at1: <Atom O>, 199, VAL
at2: <Atom N>, 203, LYS

Contact: 29
at1: <Atom O>, 206, LYS
at2: <Atom N>, 212, ALA

Contact: 30
at1: <Atom O>, 215, ALA
at2: <Atom N>, 217, LYS

Contact: 31
at1: <Atom O>, 220, LYS
at2: <Atom N>, 226, ILE

Contact: 32
at1: <Atom O>, 229, ILE
at2: <Atom N>, 234, GLN

Contact: 33
at1: <Atom O>, 237, GLN
at2: <Atom N>, 243, ASP

Contact: 34
at1: <Atom O>, 246, ASP
at2: <Atom N>, 251, LYS

Contact: 35
at1: <Atom O>, 254, LYS
at2: <Atom N>, 260, GLU

Contact: 36
at1: <Atom O>, 263, GLU
at2: <Atom N>, 269, GLY

Contact: 37
at1: <Atom O>, 272, GLY
at2: <Atom N>, 273, ILE

Contact: 38
at1: <Atom O>, 272, GLY
at2: <Atom O>, 632, HOH

Contact: 39
at1: <Atom O>, 276, ILE
at2: <Atom N>, 281, PRO

Contact: 40
at1: <Atom O>, 284, PRO
at2: <Atom N>, 288, PRO

Contact: 41
at1: <Atom O>, 291, PRO
at2: <Atom N>, 295, ASP

Contact: 42
at1: <Atom O>, 298, ASP
at2: <Atom N>, 303, GLN

Contact: 43
at1: <Atom O>, 306, GLN
at2: <Atom N>, 312, GLN

Contact: 44
at1: <Atom O>, 315, GLN
at2: <Atom N>, 321, ARG

Contact: 45
at1: <Atom O>, 324, ARG
at2: <Atom N>, 332, LEU

Contact: 46
at1: <Atom O>, 335, LEU
at2: <Atom N>, 340, ILE

Contact: 47
at1: <Atom O>, 343, ILE
at2: <Atom N>, 348, PHE

Contact: 48
at1: <Atom O>, 351, PHE
at2: <Atom N>, 359, ALA

Contact: 49
at1: <Atom O>, 362, ALA
at2: <Atom N>, 364, GLY

Contact: 50
at1: <Atom O>, 367, GLY
at2: <Atom N>, 368, LYS

Contact: 51
at1: <Atom O>, 371, LYS
at2: <Atom N>, 377, GLN

Contact: 52
at1: <Atom O>, 380, GLN
at2: <Atom N>, 386, LEU

Contact: 53
at1: <Atom O>, 389, LEU
at2: <Atom N>, 394, GLU

Contact: 54
at1: <Atom O>, 397, GLU
at2: <Atom N>, 403, ASP

Contact: 55
at1: <Atom O>, 406, ASP
at2: <Atom N>, 411, GLY

Contact: 56
at1: <Atom O>, 414, GLY
at2: <Atom N>, 415, ARG

Contact: 57
at1: <Atom O>, 418, ARG
at2: <Atom N>, 426, THR

Contact: 58
at1: <Atom O>, 429, THR
at2: <Atom N>, 433, LEU

Contact: 59
at1: <Atom O>, 436, LEU
at2: <Atom N>, 441, SER

Contact: 60
at1: <Atom O>, 444, SER
at2: <Atom N>, 447, ASP

Contact: 61
at1: <Atom O>, 450, ASP
at2: <Atom N>, 455, TYR

Contact: 62
at1: <Atom O>, 458, TYR
at2: <Atom N>, 467, ASN

Contact: 63
at1: <Atom O>, 470, ASN
at2: <Atom N>, 475, ILE

Contact: 64
at1: <Atom O>, 478, ILE
at2: <Atom N>, 483, GLN

Contact: 65
at1: <Atom O>, 486, GLN
at2: <Atom N>, 492, LYS

Contact: 66
at1: <Atom O>, 495, LYS
at2: <Atom N>, 501, GLU

Contact: 67
at1: <Atom O>, 504, GLU
at2: <Atom N>, 510, SER

Contact: 68
at1: <Atom O>, 513, SER
at2: <Atom N>, 516, THR

Contact: 69
at1: <Atom O>, 519, THR
at2: <Atom N>, 523, LEU

Contact: 70
at1: <Atom O>, 526, LEU
at2: <Atom N>, 531, HIS

Contact: 71
at1: <Atom O>, 534, HIS
at2: <Atom N>, 541, LEU

Contact: 72
at1: <Atom O>, 544, LEU
at2: <Atom N>, 549, VAL

Contact: 73
at1: <Atom O>, 552, VAL
at2: <Atom N>, 556, LEU

Contact: 74
at1: <Atom O>, 559, LEU
at2: <Atom N>, 564, ARG

Contact: 75
at1: <Atom O>, 567, ARG
at2: <Atom N>, 575, LEU

Contact: 76
at1: <Atom O>, 578, LEU
at2: <Atom N>, 583, ARG

Contact: 77
at1: <Atom O>, 586, ARG
at2: <Atom N>, 594, GLY

Contact: 78
at1: <Atom O>, 597, GLY
at2: <Atom N>, 598, GLY

Contact: 79
at1: <Atom O>, 641, HOH
at2: <Atom O>, 649, HOH

Contact: 80
at1: <Atom O>, 648, HOH
at2: <Atom O>, 655, HOH

3) Invalid input - negative distance: python 'argparse(3).py' 1ubq.pdb -d -1

Settings
--------
distance  : -1.0
PDB_file  : 1ubq.pdb

Settings, again
---------------
-1.0 1ubq.pdb
ATOM: MET, 1, N
ATOM: MET, 1, O
ATOM: MET, 1, SD
ATOM: GLN, 2, N
ATOM: GLN, 2, O
ATOM: GLN, 2, OE1
ATOM: GLN, 2, NE2
ATOM: ILE, 3, N
ATOM: ILE, 3, O
ATOM: PHE, 4, N
ATOM: PHE, 4, O
ATOM: VAL, 5, N
ATOM: VAL, 5, O
ATOM: LYS, 6, N
ATOM: LYS, 6, O
ATOM: LYS, 6, NZ
ATOM: THR, 7, N
ATOM: THR, 7, O
ATOM: THR, 7, OG1
ATOM: LEU, 8, N
ATOM: LEU, 8, O
ATOM: THR, 9, N
ATOM: THR, 9, O
ATOM: THR, 9, OG1
ATOM: GLY, 10, N
ATOM: GLY, 10, O
ATOM: LYS, 11, N
ATOM: LYS, 11, O
ATOM: LYS, 11, NZ
ATOM: THR, 12, N
ATOM: THR, 12, O
ATOM: THR, 12, OG1
ATOM: ILE, 13, N
ATOM: ILE, 13, O
ATOM: THR, 14, N
ATOM: THR, 14, O
ATOM: THR, 14, OG1
ATOM: LEU, 15, N
ATOM: LEU, 15, O
ATOM: GLU, 16, N
ATOM: GLU, 16, O
ATOM: GLU, 16, OE1
ATOM: GLU, 16, OE2
ATOM: VAL, 17, N
ATOM: VAL, 17, O
ATOM: GLU, 18, N
ATOM: GLU, 18, O
ATOM: GLU, 18, OE1
ATOM: GLU, 18, OE2
ATOM: PRO, 19, N
ATOM: PRO, 19, O
ATOM: SER, 20, N
ATOM: SER, 20, O
ATOM: SER, 20, OG
ATOM: ASP, 21, N
ATOM: ASP, 21, O
ATOM: ASP, 21, OD1
ATOM: ASP, 21, OD2
ATOM: THR, 22, N
ATOM: THR, 22, O
ATOM: THR, 22, OG1
ATOM: ILE, 23, N
ATOM: ILE, 23, O
ATOM: GLU, 24, N
ATOM: GLU, 24, O
ATOM: GLU, 24, OE1
ATOM: GLU, 24, OE2
ATOM: ASN, 25, N
ATOM: ASN, 25, O
ATOM: ASN, 25, OD1
ATOM: ASN, 25, ND2
ATOM: VAL, 26, N
ATOM: VAL, 26, O
ATOM: LYS, 27, N
ATOM: LYS, 27, O
ATOM: LYS, 27, NZ
ATOM: ALA, 28, N
ATOM: ALA, 28, O
ATOM: LYS, 29, N
ATOM: LYS, 29, O
ATOM: LYS, 29, NZ
ATOM: ILE, 30, N
ATOM: ILE, 30, O
ATOM: GLN, 31, N
ATOM: GLN, 31, O
ATOM: GLN, 31, OE1
ATOM: GLN, 31, NE2
ATOM: ASP, 32, N
ATOM: ASP, 32, O
ATOM: ASP, 32, OD1
ATOM: ASP, 32, OD2
ATOM: LYS, 33, N
ATOM: LYS, 33, O
ATOM: LYS, 33, NZ
ATOM: GLU, 34, N
ATOM: GLU, 34, O
ATOM: GLU, 34, OE1
ATOM: GLU, 34, OE2
ATOM: GLY, 35, N
ATOM: GLY, 35, O
ATOM: ILE, 36, N
ATOM: ILE, 36, O
ATOM: PRO, 37, N
ATOM: PRO, 37, O
ATOM: PRO, 38, N
ATOM: PRO, 38, O
ATOM: ASP, 39, N
ATOM: ASP, 39, O
ATOM: ASP, 39, OD1
ATOM: ASP, 39, OD2
ATOM: GLN, 40, N
ATOM: GLN, 40, O
ATOM: GLN, 40, OE1
ATOM: GLN, 40, NE2
ATOM: GLN, 41, N
ATOM: GLN, 41, O
ATOM: GLN, 41, OE1
ATOM: GLN, 41, NE2
ATOM: ARG, 42, N
ATOM: ARG, 42, O
ATOM: ARG, 42, NE
ATOM: ARG, 42, NH1
ATOM: ARG, 42, NH2
ATOM: LEU, 43, N
ATOM: LEU, 43, O
ATOM: ILE, 44, N
ATOM: ILE, 44, O
ATOM: PHE, 45, N
ATOM: PHE, 45, O
ATOM: ALA, 46, N
ATOM: ALA, 46, O
ATOM: GLY, 47, N
ATOM: GLY, 47, O
ATOM: LYS, 48, N
ATOM: LYS, 48, O
ATOM: LYS, 48, NZ
ATOM: GLN, 49, N
ATOM: GLN, 49, O
ATOM: GLN, 49, OE1
ATOM: GLN, 49, NE2
ATOM: LEU, 50, N
ATOM: LEU, 50, O
ATOM: GLU, 51, N
ATOM: GLU, 51, O
ATOM: GLU, 51, OE1
ATOM: GLU, 51, OE2
ATOM: ASP, 52, N
ATOM: ASP, 52, O
ATOM: ASP, 52, OD1
ATOM: ASP, 52, OD2
ATOM: GLY, 53, N
ATOM: GLY, 53, O
ATOM: ARG, 54, N
ATOM: ARG, 54, O
ATOM: ARG, 54, NE
ATOM: ARG, 54, NH1
ATOM: ARG, 54, NH2
ATOM: THR, 55, N
ATOM: THR, 55, O
ATOM: THR, 55, OG1
ATOM: LEU, 56, N
ATOM: LEU, 56, O
ATOM: SER, 57, N
ATOM: SER, 57, O
ATOM: SER, 57, OG
ATOM: ASP, 58, N
ATOM: ASP, 58, O
ATOM: ASP, 58, OD1
ATOM: ASP, 58, OD2
ATOM: TYR, 59, N
ATOM: TYR, 59, O
ATOM: TYR, 59, OH
ATOM: ASN, 60, N
ATOM: ASN, 60, O
ATOM: ASN, 60, OD1
ATOM: ASN, 60, ND2
ATOM: ILE, 61, N
ATOM: ILE, 61, O
ATOM: GLN, 62, N
ATOM: GLN, 62, O
ATOM: GLN, 62, OE1
ATOM: GLN, 62, NE2
ATOM: LYS, 63, N
ATOM: LYS, 63, O
ATOM: LYS, 63, NZ
ATOM: GLU, 64, N
ATOM: GLU, 64, O
ATOM: GLU, 64, OE1
ATOM: GLU, 64, OE2
ATOM: SER, 65, N
ATOM: SER, 65, O
ATOM: SER, 65, OG
ATOM: THR, 66, N
ATOM: THR, 66, O
ATOM: THR, 66, OG1
ATOM: LEU, 67, N
ATOM: LEU, 67, O
ATOM: HIS, 68, N
ATOM: HIS, 68, O
ATOM: HIS, 68, ND1
ATOM: HIS, 68, NE2
ATOM: LEU, 69, N
ATOM: LEU, 69, O
ATOM: VAL, 70, N
ATOM: VAL, 70, O
ATOM: LEU, 71, N
ATOM: LEU, 71, O
ATOM: ARG, 72, N
ATOM: ARG, 72, O
ATOM: ARG, 72, NE
ATOM: ARG, 72, NH1
ATOM: ARG, 72, NH2
ATOM: LEU, 73, N
ATOM: LEU, 73, O
ATOM: ARG, 74, N
ATOM: ARG, 74, O
ATOM: ARG, 74, NE
ATOM: ARG, 74, NH1
ATOM: ARG, 74, NH2
ATOM: GLY, 75, N
ATOM: GLY, 75, O
ATOM: GLY, 76, N
ATOM: GLY, 76, O
ATOM: GLY, 76, OXT
ATOM: HOH, 77, O
ATOM: HOH, 78, O
ATOM: HOH, 79, O
ATOM: HOH, 80, O
ATOM: HOH, 81, O
ATOM: HOH, 82, O
ATOM: HOH, 83, O
ATOM: HOH, 84, O
ATOM: HOH, 85, O
ATOM: HOH, 86, O
ATOM: HOH, 87, O
ATOM: HOH, 88, O
ATOM: HOH, 89, O
ATOM: HOH, 90, O
ATOM: HOH, 91, O
ATOM: HOH, 92, O
ATOM: HOH, 93, O
ATOM: HOH, 94, O
ATOM: HOH, 95, O
ATOM: HOH, 96, O
ATOM: HOH, 97, O
ATOM: HOH, 98, O
ATOM: HOH, 99, O
ATOM: HOH, 100, O
ATOM: HOH, 101, O
ATOM: HOH, 102, O
ATOM: HOH, 103, O
ATOM: HOH, 104, O
ATOM: HOH, 105, O
ATOM: HOH, 106, O
ATOM: HOH, 107, O
ATOM: HOH, 108, O
ATOM: HOH, 109, O
ATOM: HOH, 110, O
ATOM: HOH, 111, O
ATOM: HOH, 112, O
ATOM: HOH, 113, O
ATOM: HOH, 114, O
ATOM: HOH, 115, O
ATOM: HOH, 116, O
ATOM: HOH, 117, O
ATOM: HOH, 118, O
ATOM: HOH, 119, O
ATOM: HOH, 120, O
ATOM: HOH, 121, O
ATOM: HOH, 122, O
ATOM: HOH, 123, O
ATOM: HOH, 124, O
ATOM: HOH, 125, O
ATOM: HOH, 126, O
ATOM: HOH, 127, O
ATOM: HOH, 128, O
ATOM: HOH, 129, O
ATOM: HOH, 130, O
ATOM: HOH, 131, O
ATOM: HOH, 132, O
ATOM: HOH, 133, O
ATOM: HOH, 134, O
NBSEARCH:
Traceback (most recent call last):
  File "/Users/darialunina/Documents/Biophysics_exercises/argparse(3).py", line 69, in <module>
    for a, b in nbsearch.search_all(MAXDIST)
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/darialunina/miniconda3/lib/python3.12/site-packages/Bio/PDB/NeighborSearch.py", line 115, in search_all
    neighbors = self.kdt.neighbor_search(radius)
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ValueError: Radius must be positive.





"""