/*=============================================================================
#     FileName: atom_type.h
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-03-22 16:32:30
#   LastChange: 2014-06-19 15:34:45
#      History:
=============================================================================*/

#ifndef  ATOM_TYPE_H
#define  ATOM_TYPE_H

#include <vector>
#include <string>
#include <openbabel/mol.h>
#include <openbabel/atom.h>

/*
the following atom types are based on:
- Muegge I. and Martin Y.C. J. Med. Chem. 1999, 42: 791-804
- Muegge I. J. Med. Chem. 2006, 49:5895-5902
- Shen Q.C. et al. J. Chem. Inf. Model. 2010, 51:386-397

1. LIGAND ATOM TYPE (31)
CF: nonpolar carbon sp3 aliphatic
CP: polar sp3 carbon bonded to an atom other than carbon or hydrogen
cF: nonpolar carbon aromatic
cP: polar carbon aromatic
C3: nonpolar carbon sp2 not aromatic
CW: polar carbon sp2 not aromatic (e.g., bonded to carbonyl oxygen)
CO: carbon bonded to a negatively charged oxygen
CN: carbon bonded to a positively charged nitrogen
C0: sp carbon
NC: positively charged nitrogen (e.g., NH3+ or guanidino group)
NP: planar nitrogen bonded to 2 or 3 carbons but not to a hydrogen
    (can occur in a nonaromatic ring)
NA: nitrogen as a hydrogen bond acceptor, not in a ring
ND: nitrogen as a hydrogen bond donor, not in a ring (e.g., amide nitrogen)
NR: planar nitrogen in a ring structure (e.g., pyridine)
N0: sp nitrogen bonded to 1 carbon
NS: nitrogen boound to atoms other than carbon or hydrogen and not of type ND
OC: negatively charged oxygen (e.g., carboxylate)
OA: oxygen as hydrogen bond acceptor (e.g., keto, amide oxygen)
OE: oxygen in an ether bond or ring
OS: oxygen bonded to atoms other than carbon or hydrogen (e.g. phosphate)
OD: oxygen bonded to hydrogen, except water
P : phosphorus
SA: sulfur as hydrogen bond acceptor (except solfune/sulfoxide)
SD: sulfur as hydrogen bond donor
SO: sulfur in sulfone or sulfoxide
HL: hydrogen
F : fluorine
CL: chlorine
BR: bromine
I : iodine
ME: metal (Zn, Mn, Mg, Fe, V)

2. PROTEIN ATOM TYPE (17)
CF: nonpolar aliphatic carbon (e.g., CB)
CP: polar aliphatic SP2 or SP3 carbon bonded to atoms other than carbon or hydrogen
    (e.g., backbone C or Ca)
cF: nonpolar carbon aromatic
cP: polar carbon aromatic
CO: carbon bonded to a negatively charged oxygen
CN: carbon bonded to a positively charged nitrogen
NA: nitrogen as hydrogen bond acceptor (e.g. HIS NE2)
NC: positively charged nitrogen
ND: nitrogen as hydrogen bond donor (e.g., backbone N, TRP NE, ASN ND, HIS ND1)
OC: negatively charged oxygen
OA: oxygen as hydrogen bond acceptor (e.g., backbone O, ASN OD, GLN OE)
OD: oxygen as hydrogen bond donor (e.g., TYR OH, SER OG, THR OG)
OW: water oxygen
SA: sulfur as hydrogen bond acceptor (MET SD)
SD: sulfur as hydrogen bond donor (CYS SG)
HH: hydrogen
ME: metal, including Zn, Ca, K, Mg, Mn, and Fe
*/
/* 
LTYPE = [
        1~9  : "CF","CP","cF","cP","C3","CW","CO","CN","C0",
        10~16: "NC","NP","NA","ND","NR","N0","NS",
        17~21: "OC","OA","OE","OS","OD",
        22~31: "P","SA","SD","SO","HL","F","CL","BR","I","ME"]
PTYPE = [
        1~6  : "CF","CP","cF","cP","CO","CN",
        7~9  : "NA","NC","ND",
        10~13: "OC","OA","OD","OW",
        14~17: "SA","SD","HH","ME"
        ]
*/
extern std::vector<std::string> LTYPE;
extern std::vector<std::string> PTYPE;

// if no pre-defined atom type is matched, return 0
int parse_ligand_atom_v1(OpenBabel::OBAtom &atom);
int* parse_ligand_v1(OpenBabel::OBMol &ligand);

// if no pre-defined atom type is matched, return 0
int parse_protein_atom_v1(OpenBabel::OBAtom &atom);
int* parse_protein_v1(OpenBabel::OBMol &protein);

#endif   /* ----- #ifndef ATOM_TYPE_H  ----- */

