'''
#=============================================================================
#     FileName: atom_type.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-03-21 09:29:50
#   LastChange: 2014-07-12 01:05:36
#      History:
#=============================================================================
'''
import sys
import pybel

pmf04_atom_type = """
the following atom types are based on:
- Muegge I. and Martin Y.C. J. Med. Chem. 1999, 42: 791-804
- Muegge I. J. Med. Chem. 2006, 49:5895-5902
- Shen Q.C. et al. J. Chem. Inf. Model. 2010, 51:386-397

1. LIGAND ATOM TYPE (30)
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
#ME: metal (Zn, Mn, Mg, Fe, V, Co)

2. PROTEIN ATOM TYPE (17)
CF: nonpolar aliphatic carbon (e.g., CB)
CP: polar aliphatic SP2 or SP3 carbon bonded to atoms other than carbon or
    hydrogen (e.g., backbone C or Ca)
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
"""

LTYPE = [("CF","CP","cF","cP","C3","CW","CO","CN","C0",\
        "NC","NP","NA","ND","NR","N0","NS",\
        "OC","OA","OE","OS","OD",\
        "P","SA","SD","SO","HL","F","CL","BR","I","ME")]
PTYPE = [("CF","CP","cF","cP","CO","CN",\
        "NA","NC","ND",\
        "OC","OA","OD","OW",\
        "SA","SD","HH","ME")]

atom_type_description = [pmf04_atom_type]

#to specify the available '_type'
#please update this tuple whenever you add more types in LTYPE and PTYPE
available_types = (0,)


def numConnectH(atom):
    num_h = 0
    for nbor in pybel.ob.OBAtomAtomIter(atom):
        if nbor.IsHydrogen():
            num_h += 1
    return num_h

def isChargedNitrogen(atom):
    atype = atom.GetType()
    isCharged = False
    if(not atom.IsNitrogen()):
        return False
    if(atype == "N3"):
        num_het = 0
        for nbor in pybel.ob.OBAtomAtomIter(atom):
            nbor_type = nbor.GetType()
            if nbor_type!="C3" and (not nbor.IsHydrogen()):
                num_het += 1
        if num_het == 0:
            isCharged = True
    elif atype == "Ng+":
        if atom.ImplicitHydrogenCount()==2 or atom.ExplicitHydrogenCount()==2:
            isCharged = True
    elif atom.KBOSum() == 4:
        isCharged = True
    else:
        isCharged = False
    return isCharged


def isChargedOxygen(atom):
    atype = atom.GetType()
    if atype in ["O-","O.co2","OCO2"] or \
            atom.MatchesSMARTS("[$([#8-]),$([OX2H1]C=O),$(O=C[OX2H1])]"):
        return True
    else:
        return False


"""
the following atom types are based on PMF99
"""
#INPUT: OBAtom
def parse_ligand_atom_v1(atom):
    ltype = ""
    atomic_num = atom.GetAtomicNum()
    if atomic_num == 6:
        num_connect_het = 0
        for nbor in pybel.ob.OBAtomAtomIter(atom):
            if (not nbor.IsCarbon()) and (not nbor.IsHydrogen()):
                num_connect_het += 1
            if isChargedOxygen(nbor):
                ltype = "CO"
                break
            if isChargedNitrogen(nbor):
                ltype = "CN"
                break
        if ltype != "":
            return ltype
        if atom.IsAromatic():
            #if atom.MatchesSMARTS("c~[!#6;!#1]"):
            if num_connect_het:
                ltype = "cP"
            else:
                ltype = "cF"
        elif atom.GetHyb() == 3:
            #if atom.MatchesSMARTS("[C^3]~[!#6;!#1]"):
            if num_connect_het:
                ltype = "CP"
            else:
                ltype = "CF"
        elif atom.GetHyb() == 2:
            #if atom.MatchesSMARTS("[C^2]~[!#6;!#1]"):
            if num_connect_het:
                ltype = "CW"
            else:
                ltype = "C3"
        elif atom.GetHyb() == 1:
            ltype = "C0"
        else:
            pass
    elif atomic_num == 7:
        num_connect_h = atom.ImplicitHydrogenCount()
        num_connect_c = 0
        num_connect_het = 0
        for nbor in pybel.ob.OBAtomAtomIter(atom):
            if nbor.IsHydrogen():
                num_connect_h += 1
            elif nbor.IsCarbon():
                num_connect_c += 1
            else:
                num_connect_het += 1
        if not atom.IsInRing():
            if atom.IsHbondAcceptor() or (num_connect_h==0):
                ltype = "NA"
            if atom.IsHbondDonor() or (num_connect_h>0):
                ltype = "ND"
        else:
            if atom.GetHyb() == 2:
                ltype = "NR"
        if num_connect_c>=2 and num_connect_h==0 and (not atom.IsAromatic()) and \
                (atom.IsInRing()) and (atom.GetHyb()==2):
            ltype = "NP"
        if num_connect_c==1 and atom.GetHyb()==1:
            ltype = "N0"    # "[#7^1]#[#6]"
        if isChargedNitrogen(atom):
            ltype = "NC"
        if num_connect_het>0 and ltype!="ND":
            ltype = "NS"
    elif atomic_num == 8:
        if isChargedOxygen(atom):
            ltype = "OC"
        else:
            num_connect_c = 0
            num_connect_het = 0
            num_connect_h = atom.ImplicitHydrogenCount()
            for nbor in pybel.ob.OBAtomAtomIter(atom):
                if nbor.IsHydrogen():
                    num_connect_h += 1
                elif nbor.IsCarbon():
                    num_connect_c += 1
                else:
                    num_connect_het += 1
            if atom.IsHbondAcceptor():
                ltype = "OA"
            #if atom.IsHbondDonor() or atom.MatchesSMARTS("[#8;!H0]"):
            if atom.IsHbondDonor() or num_connect_h>0:
                ltype = "OD"
            #if atom.MatchesSMARTS("[#8](~[#6])~[#6]"):
            if num_connect_c==2 or atom.IsInRing():
                ltype = "OE"
            #if atom.MatchesSMARTS("[#8]~[!#6;!#1]"):
            if num_connect_het:
                ltype = "OS"
    elif atomic_num == 15:
        ltype = "P"
    elif atomic_num == 16:
        ltype = "SA"
        #if atom.IsHbondDonor() or atom.MatchesSMARTS("[#16;!H0]"):
        if atom.IsHbondDonor() or atom.ImplicitHydrogenCount()>0 or atom.ExplicitHydrogenCount()>0:
            ltype = "SD"
        if atom.MatchesSMARTS("[#16]=[#8]"):
            ltype = "SO"
    elif atomic_num == 1:
        ltype = "HL";
    elif atomic_num == 9:
        ltype = "F"
    elif atomic_num == 17:
        ltype = "CL"
    elif atomic_num == 35:
        ltype = "BR"
    elif atomic_num == 53:
        ltype = "I"
    elif atomic_num in [30,25,12,26,23,27]:    #Zn, Mn, Mg, Fe, V, Co
        ltype = "ME"
    else:
        ltype = ""

    return ltype


def parse_protein_atom_v1(atom):
    ptype = ""
    atomic_num = atom.GetAtomicNum()
    if atomic_num == 6:
        num_connect_het = 0
        for nbor in pybel.ob.OBAtomAtomIter(atom):
            if((not nbor.IsCarbon()) and (not nbor.IsHydrogen())):
                num_connect_het += 1
            if isChargedOxygen(nbor):
                ptype = "CO"
                break
            if isChargedNitrogen(nbor):
                ptype = "CN"
                break
        if ptype != "":
            return ptype
        if atom.IsAromatic():
            if num_connect_het:
            #if atom.MatchesSMARTS("c~[!#6;!#1]"):
                ptype = "cP"
            else:
                ptype = "cF"
        else:
            if num_connect_het:
            #if atom.MatchesSMARTS("[#6]~[!#6;!#1]"):
                ptype = "CP"
            else:
                ptype = "CF"
    elif atomic_num == 7:
        ptype = "NA"
        #if atom.IsHbondAcceptor():
        #    ptype = "NA"
        #if atom.IsHbondDonor() or atom.MatchesSMARTS("[#7;!H0]"):
        if atom.IsHbondDonor() or atom.ImplicitHydrogenCount()>0 or atom.ExplicitHydrogenCount()>0:
            ptype = "ND"
        if isChargedNitrogen(atom):
            ptype = "NC"
    elif atomic_num == 8:
        if atom.MatchesSMARTS("[#8;H2]"):
            ptype = "OW"
            return ptype
        ptype = "OA"
        #if atom.IsHbondAcceptor():
        #    ptype = "OA"
        #if atom.IsHbondDonor() or atom.MatchesSMARTS("[#8;!H0]"):
        if atom.IsHbondDonor() or atom.ImplicitHydrogenCount()>0 or atom.ExplicitHydrogenCount()>0:
            ptype = "OD"
        if isChargedOxygen(atom):
            ptype = "OC"
    elif atomic_num == 16:
        ptype = "SA"
        #if atom.IsHbondDonor() or atom.MatchesSMARTS("[#16;!H0]"):
        if atom.IsHbondDonor() or atom.ImplicitHydrogenCount()>0 or atom.ExplicitHydrogenCount()>0:
            ptype = "SD"
    elif atomic_num == 1:
        ptype = "HH"
    elif atomic_num in [30,12,20,26,25,19]:   #Zn, Mg, Ca, Fe, Mn, K
        ptype = "ME"
    else:
        ptype = ""
    
    return ptype



"""
the following atom types are based on M-Score:
- 
"""
def parse_ligand_atom_v2(atom):
    ltype = ""
    atomic_num = atom.GetAtomicNum()
    hyb = atom.GetHyb()
    if atomic_num == 6:
        if atom.MatchesSMARTS("[#6+]"):
            ltype = "C.cat"
        elif hyb == 3:
            ltype = "C.3"
        elif hyb == 2:
            if atom.IsAromatic():
                ltype = "C.ar"
            else:
                ltype = "C.2"
        elif hyb == 1:
            ltype = "C.1"
        else:
            pass
    elif atomic_num == 7:
        #N.pl3  ???
        if atom.IsAromatic():
            ltype = "N.ar"
        elif atom.IsAmideNitrogen():
            ltype = "N.am"
        elif hyb == 3:
            if atom.MatchesSMARTS("[#7+]"):
                ltype = "N.4"
            else:
                ltype = "N.3"
        elif hyb == 2:
            ltype = "N.2"
        elif hyb == 1:
            ltype = "N.1"
        else:
            pass
    elif atomic_num == 8:
        if hyb == 3:
            ltype = "O.3"
        elif hyb == 2:
            if atom.IsCarboxylOxygen() or atom.IsPhosphateOxygen():
                ltype = "O.co2"
            else:
                ltype = "O.2"
        else:
            pass
    elif atomic_num == 15:
        ltype = "P.3"
    elif atomic_num == 16:
        if hyb == 3:
            ltype = "S.3"
        elif hyb == 2:
            if atom.MatchesSMARTS("[#16]=O"):
                ltype = "S.O"
            else:
                ltype = "S.2"
        else:
            pass
    elif atomic_num == 1:
        ltype = "H"
    elif atomic_num == 9:
        ltype = "F"
    elif atomic_num == 17:
        ltype = "CL"
    elif atomic_num == 35:
        ltype = "BR"
    elif atomic_num == 53:
        ltype = "I"
    else:
        ltype = "ME"
    
    return ltype
    

def parse_protein_atom_v2(atom):
    ptype = ""
    atomic_num = atom.GetAtomicNum()
    res = a.GetResidue().GetName()
    if atomic_num == 6:
        if atom.MatchesSMARTS("C(N)C=O"):  #alpha carbon
            ptype = "Ca"
        elif atom.MatchesSMARTS("C(=O)CN"): #backbone carbonyl carbon
            ptype = "Cbk"
        elif atom.MatchesSMARTS("CC(N)C=O"):  #beta carbon
            ptype = "Cb"
        elif atom.GetHyb == 3:
            if atom.MatchesSMARTS("C(~N)(~N)~N"):  #charged carbon in R
                ptype = "C.cat(R)"
            else:                                  #side chain sp3 carbon
                ptype = "C.3(sch)"
        elif atom.GetHyb == 2:
            if atom.IsAromatic():  #carbon in aromatic system
                ptype = "C.ar"
            else:                  #side chain sp2 carbon
                if res in ["ASN","GLN"]:
                    ptype = "C.2(N,Q)"
                elif res in ["ASP","GLU"]:
                    ptype = "C.2(E,D)"
                else:
                    pass
        else:
            pass
    elif atomic_num == 7:
        if atom.MatchesSMARTS("NCC=O"):  #backbone amide nitrogen
            ptype = "N.am"
        elif atom.IsAromatic():          #nitrogen in aromatic system (W, H)
            ptype = "N.ar"
        elif res in ["ASN","GLN"]:       #polar and uncharged nitrogen (N, Q)
            ptype = "Npolar(N,Q)"
        elif res in ["ARG","LYS"]:       #charged nitrogen atom (R, K)
            ptype = "N.4(R,K)"
        else:
            pass
    elif atomic_num == 8:
        if atom.MatchesSMARTS("O=CCN"):  #backbone carbonyl oxygen
            ptype = "O(bk)"
        elif atom.GetHyb() == 2:
            if res in ["GLU","ASP"]:     #carboxylic charged oxygen (E,D)
                ptype = "O.co2"
            elif res in ["GLN","ASN"]:   #carboxylic uncharged oxygen (N,Q)
                ptype = "O.2(N,Q)"
            else:
                pass
        elif atom.MatchesSMARTS("[OX2H1]"): #polar and not charged oxygen
            ptype = "O(polar)"
        else:
            pass
    elif atomic_num == 16:
        if res == "MET":
            ptype = "S(M)"
        elif res == "CYS":
            ptype = "S(C)"
        else:
            pass
    elif atomic_num == 1:
        ptype = "HH"
    else:
        ptype = "ME"

    return ptype


#=================================================
#please use 'parse_ligand' and 'parse_protein' instead
#the order must be compatible with 'PTYPE' and 'LTYPE'
parse_ligand_atom = [parse_ligand_atom_v1, parse_ligand_atom_v2]
parse_protein_atom = [parse_protein_atom_v1, parse_protein_atom_v2]

#INPUT:
#  ligand - OBMol
#  _type  - int, which function of 'parse_ligand_atom' to be used
def parse_ligand(ligand, _type):
    ltypes = ["" for _ in xrange(ligand.NumAtoms())]
    for atom in pybel.ob.OBMolAtomIter(ligand):
        ltypes[atom.GetIdx()-1] = parse_ligand_atom[_type](atom)
    return ltypes

#parse_ligand_v1 = lambda ligand: [parse_ligand_atom[_type](atom) for atom in pybel.ob.OBMolAtomIter(ligand)]

#INPUT:
#  protein - OBMol
#  _type   - int, which function of 'parse_protein_atom' to be used
def parse_protein(protein, _type):
    ptypes = ["" for _ in xrange(protein.NumAtoms())]
    for atom in pybel.ob.OBMolAtomIter(protein):
        ptypes[atom.GetIdx()-1] = parse_protein_atom[_type](atom)
    return ptypes

#parse_protein_v1 = lambda protein: [parse_protein_atom[_type](atom) for atom in pybel.ob.OBMolAtomIter(protein)]
