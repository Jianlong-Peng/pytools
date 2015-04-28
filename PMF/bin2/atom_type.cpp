/*=============================================================================
#     FileName: atom_type.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-03-22 16:32:48
#   LastChange: 2014-08-15 11:48:11
#      History:
=============================================================================*/
#include <vector>
#include <string>
#include <cstring>
#include "atom_type.h"
#include <openbabel/mol.h>
#include <openbabel/atom.h>

using namespace std;
using namespace OpenBabel;

string temp_ltype[] = {
        /* 1~9   */"CF","CP","AF","AP","C3","CW","CO","CN","C0",
        /* 10~16 */"NC","NP","NA","ND","NR","N0","NS",
        /* 17~21 */"OC","OA","OE","OS","OD",
        /* 22~30 */"P","SA","SD","SO","HL","F","CL","Br","I"};
vector<string> LTYPE(temp_ltype, temp_ltype+30);
string temp_ptype[] = {
    /* 1~6   */"CF","CP","AF","AP","CO","CN",
    /* 7~9   */"NA","NC","ND",
    /* 10~13 */"OC","OA","OD","OW",
    /* 14~17 */"SA","SD","HH","ME"};
vector<string> PTYPE(temp_ptype, temp_ptype+17);


inline static int numConnectH(OBAtom &atom)
{
    int num_connected_h = 0;
    FOR_NBORS_OF_ATOM(nbor, atom) {
        if(nbor->IsHydrogen())
            ++num_connected_h;
    }
    return num_connected_h;
}

/*
inline static bool isChargedNitrogen(OBAtom &atom)
{
    if(atom.MatchesSMARTS("[$([#7;+]),$([#7X4])]") || strcmp(atom.GetType(),"Ng+")==0)
        return true;
    else
        return false;
}
*/
static bool isChargedNitrogen(OBAtom &atom)
{
    char *type = atom.GetType();
    bool isCharged = false;
    if(!atom.IsNitrogen())
        return false;
    if(strcmp(type,"N3")==0) {   // nitrogen is supposed to be protonated!!!
        int num_het = 0;
        FOR_NBORS_OF_ATOM(nbor, atom) {
            char *nbor_type = nbor->GetType();
            if(strcmp(nbor_type,"C3") && !nbor->IsHydrogen())
                ++num_het;
        }
        if(num_het == 0)
            isCharged = true;
    }
    else if(strcmp(type, "Ng+")==0) {
        if(atom.ImplicitHydrogenCount()==2 || atom.ExplicitHydrogenCount()==2)
            isCharged = true;
    }
    else if(atom.KBOSum() == 4)
        isCharged = true;
    else
        isCharged = false;
    
    return isCharged;
}
/*
inline static bool isChargedOxygen(OBAtom &atom)
{
    if(atom.MatchesSMARTS("[$([#8;-]),$([OX2H1]C=O),$(O=C[OX2H1])]"))
        return true;
    else
        return false;
}
*/
static bool isChargedOxygen(OBAtom &atom)
{
    char *type = atom.GetType();
    bool isCharged = false;
    if(strcmp(type,"O-")==0 || strcmp(type,"O.co2")==0 || strcmp(type,"OCO2")==0 ||
        atom.MatchesSMARTS("[$([#8-]),$([OX2H1]C=O),$(O=C[OX2H1])]"))
        isCharged = true;
    return isCharged;
}

int parse_ligand_atom_v1(OBAtom &atom)
{
    int type = 0;
    unsigned atomic_num = atom.GetAtomicNum();
    unsigned num_connected_het = atom.GetHeteroValence();
    if(atomic_num == 1)
        type = 26;
    else if(atomic_num == 6) {
        if(atom.IsAromatic()) {
            if(num_connected_het > 0)  // AP
                type = 4;
            else                       // AF
                type = 3;
        } else {
            if(atom.GetHyb() == 2) {
                if(num_connected_het > 0)  // CW
                    type = 6;
                else                       // C3
                    type = 5;
            } else if(atom.GetHyb() == 1)    // C0
                type = 9;
            else {
                if(num_connected_het > 0)  // CP
                    type = 2;
                else                       // CF
                    type = 1;
            }
        }
        FOR_NBORS_OF_ATOM(nbor, atom) {
            if(isChargedNitrogen(*nbor)) {  // CN
                type = 8;
                break;
            }
        }
        FOR_NBORS_OF_ATOM(nbor, atom) {
            if(isChargedOxygen(*nbor)) {  // CO
                type = 7;
                break;
            }
        }
    }
    else if(atomic_num == 7) {
        type = 12;  // NA
        unsigned num_connected_h = atom.ImplicitHydrogenCount();
        unsigned num_connected_c = 0;
        unsigned num_connected_het = 0;
        FOR_NBORS_OF_ATOM(nbor,atom) {
            if(nbor->IsHydrogen())
                num_connected_h += 1;
            else if(nbor->IsCarbon())
                num_connected_c += 1;
            else
                num_connected_het += 1;
        }
        if(num_connected_c==1 && atom.GetHyb()==1)  // N0
            type = 15;
        if(num_connected_het > 0)  // NS
            type = 16;
        if(!atom.IsInRing()) {
            if(atom.IsHbondAcceptor() || num_connected_h == 0)  // NA
                type = 12;
            if(atom.IsHbondDonor() || num_connected_h>0)  // ND
                type = 13;
        } else {
            if(atom.GetHyb() == 2) {
                type = 14;   // NR
                if(num_connected_c>=2 && num_connected_h==0 && (!atom.IsAromatic())) // NP
                    type = 11;
            }
        }
        if(isChargedNitrogen(atom))  // NC
            type = 10;
    }
    else if(atomic_num == 8) {
        type = 18;  // OA
        int num_connect_c = 0;
        int num_connect_h = atom.ImplicitHydrogenCount();
        FOR_NBORS_OF_ATOM(nbor, atom) {
            if(nbor->IsCarbon())
                ++num_connect_c;
            else if(nbor->IsHydrogen())
                ++num_connect_h;
        }
        if(atom.IsInRing())  // OE
            type = 19;
        else {
            if(num_connect_c == 2)    // OE
                type = 19;
            FOR_NBORS_OF_ATOM(nbor, atom) {
                if((!nbor->IsCarbon()) and (!nbor->IsHydrogen())) {  // OS
                    type = 20;
                    break;
                }
            }
        }
        if(atom.IsHbondDonor() || num_connect_h>0)  // OD
            type = 21;
        if(isChargedOxygen(atom))   // OC
            type = 17;
    }
    else if(atomic_num == 9)
        type = 27;
    else if(atomic_num == 15)
        type = 22;
    else if(atomic_num == 16) {
        type = 23;  // SA
        int num_connected_o = 0;
        int num_connected_h = 0;
        FOR_NBORS_OF_ATOM(nbor, atom) {
            if(nbor->IsOxygen())
                ++num_connected_o;
            if(nbor->IsHydrogen())
                ++num_connected_h;
        }
        if(num_connected_o > 0)  // SO
            type = 25;
        if(atom.IsHbondDonor() || atom.ImplicitHydrogenCount()>0 || num_connected_h>0) // SD
            type = 24;
    }
    else if(atomic_num == 17)    // Cl
        type = 28;
    else if(atomic_num == 35)    // Br
        type = 29;
    else if(atomic_num == 53)    // I
        type = 30;
    else
        type = 0;

    return type;
}

int* parse_ligand_v1(OBMol &ligand)
{
    int *types = (int*)malloc(sizeof(int)*(ligand.NumAtoms()));
    for(unsigned i=0; i<ligand.NumAtoms(); ++i)
        types[i] = parse_ligand_atom_v1(*(ligand.GetAtom(i+1)));
    return types;
}

int parse_protein_atom_v1(OBAtom &atom)
{
    unsigned atomic_num = atom.GetAtomicNum();
    int type = 0;
    if(atomic_num == 6) {
        int num_connected_het = 0;
        FOR_NBORS_OF_ATOM(nbor,atom)
            if(!nbor->IsCarbon() && !nbor->IsHydrogen())
                ++num_connected_het;
        if(!atom.IsAromatic()) {
            if(num_connected_het == 0)  // CF
                type = 1;
            else                        // CP
                type = 2;
        } else {
            if(num_connected_het == 0)  // AF
                type = 3;
            else                        // AP
                type = 4;
        }
        FOR_NBORS_OF_ATOM(nbor, atom) {
            if(isChargedNitrogen(*nbor)) {  // CN
                type = 6;
                break;
            }
        }
        FOR_NBORS_OF_ATOM(nbor, atom) {
            if(isChargedOxygen(*nbor)) {  // CO
                type = 5;
                break;
            }
        }
    }
    else if(atomic_num == 7) {
        type = 7;   // NA
        int num_connected_h = atom.ImplicitHydrogenCount();
        FOR_NBORS_OF_ATOM(nbor, atom)
            if(nbor->IsHydrogen())
                ++num_connected_h;
        if(atom.IsHbondDonor() || num_connected_h>0)  // ND
            type = 9;
        if(isChargedNitrogen(atom))  // NC
            type = 8;
    }
    else if(atomic_num == 8) {
        type = 11;  // OA
        int num_connected_h = atom.ImplicitHydrogenCount();
        FOR_NBORS_OF_ATOM(nbor, atom)
            if(nbor->IsHydrogen())
                ++num_connected_h;
        if(atom.IsHbondDonor() || num_connected_h>0)  // OD
            type = 12;
        if(isChargedOxygen(atom))  // OC
            type = 10;
        string res = atom.GetResidue()->GetName();
        if(res=="HOH" || res=="WAT")  // OW
            type = 13;
    }
    else if(atomic_num == 16) {
        type = 14;    // SA
        int num_connected_h = 0;
        FOR_NBORS_OF_ATOM(nbor, atom)
            if(nbor->IsHydrogen())
                ++num_connected_h;
        if(atom.IsHbondDonor() || num_connected_h>0)  // SD
            type = 15;
    }
    else if(atomic_num == 1)   // HH
        type = 16;
    else if(atomic_num == 30)   // Zn
        type = 17;
    else if(atomic_num == 12)   // Mg
        type = 17;
    else if(atomic_num == 20)   // Ca
        type = 17;
    else if(atomic_num == 26)   // Fe
        type = 17;
    else if(atomic_num == 25)   // Mn
        type = 17;
    else if(atomic_num == 19)   // K
        type = 17;
    else
        type = 0;

    return type;
}
int* parse_protein_v1(OBMol &protein)
{
    int *types = (int*)malloc(sizeof(int)*(protein.NumAtoms()));
    FOR_ATOMS_OF_MOL(atom, protein)
        types[atom->GetIdx()-1] = parse_protein_atom_v1(*atom);
    return types;
}

