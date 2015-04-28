/*=============================================================================
#     FileName: atom_type.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-03-22 16:32:48
#   LastChange: 2014-06-20 11:06:34
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
        /* 1~9   */"CF","CP","cF","cP","C3","CW","CO","CN","C0",
        /* 10~16 */"NC","NP","NA","ND","NR","N0","NS",
        /* 17~21 */"OC","OA","OE","OS","OD",
        /* 22~31 */"P","SA","SD","SO","HL","F","CL","BR","I","ME"};
vector<string> LTYPE(temp_ltype, temp_ltype+31);
string temp_ptype[] = {
    /* 1~6   */"CF","CP","cF","cP","CO","CN",
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
    if(atomic_num == 1)
        type = 26;
    else if(atomic_num == 6) {
        int num_connect_het = 0;
        FOR_NBORS_OF_ATOM(nbor, atom) {
            if(!nbor->IsCarbon() && !nbor->IsHydrogen())
                ++num_connect_het;
            if(isChargedOxygen(*nbor)) {  // CO
                type = 7;
                break;
            }
            if(isChargedNitrogen(*nbor)) { // CN
                type = 8;
                break;
            }
        }
        if(type == 0) {
            if(atom.IsAromatic()) {
                //if(atom.MatchesSMARTS("c~[!#6;!#1]"))   // cP
                if(num_connect_het)
                    type = 4;
                else                                    // cF
                    type = 3;
            }
            else if(atom.GetHyb() == 3) {
                //if(atom.MatchesSMARTS("[C^3]~[!#6;!#1]"))  // CP
                if(num_connect_het)
                    type = 2;
                else                                       // CF
                    type = 1;
            }
            else if(atom.GetHyb() == 2) {
                //if(atom.MatchesSMARTS("[C^2]~[!#6;!#1]"))  // CW
                if(num_connect_het)
                    type = 6;
                else                                       // C3
                    type = 5;
            }
            else                                            // C0
                type = 9;
        }
    }
    else if(atomic_num == 7) {
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
        if(atom.IsInRing()) {
            if(atom.GetHyb() == 2)   // NR   ???
                type = 14;
        } else {
            if(atom.IsHbondAcceptor() || (num_connected_h == 0))  // NA
                type = 12;
            if(atom.IsHbondDonor() || (num_connected_h > 0))      // ND
                type = 13;
        }
        if(num_connected_c>=2 && num_connected_h==0 && (!atom.IsAromatic()) &&
                atom.IsInRing() && atom.GetHyb()==2)              // NP
            type = 11;
        if(num_connected_c==1 && atom.GetHyb()==1)                // N0
            type = 15;
        if(isChargedNitrogen(atom))                               // NC
            type = 10;
        if(num_connected_het>0 && type!=13)                       // NS
            type = 16;
    }
    else if(atomic_num == 8) {
        if(isChargedOxygen(atom))   // OC
            type = 17;
        else {
            int num_connect_c = 0;
            int num_connect_het = 0;
            int num_connect_h = atom.ImplicitHydrogenCount();
            FOR_NBORS_OF_ATOM(nbor, atom) {
                if(nbor->IsHydrogen())
                    ++num_connect_h;
                else if(nbor->IsCarbon())
                    ++num_connect_c;
                else
                    ++num_connect_het;
            }
            if(atom.IsHbondAcceptor())  // OA
                type = 18;
            //if(atom.IsHbondDonor() || atom.MatchesSMARTS("[#8;!H0]"))  // OD
            if(atom.IsHbondDonor() || num_connect_h>0)
                type = 21;
            //if(atom.MatchesSMARTS("[#8](~[#6])~[#6]") || atom.IsInRing())  // OE
            if(num_connect_c==2 || atom.IsInRing())
                type = 19;
            //if(atom.MatchesSMARTS("[#8]~[!#6;!#1]"))    // OS
            if(num_connect_het)
                type = 20;
        }
    }
    else if(atomic_num == 9)
        type = 27;
    else if(atomic_num == 15)
        type = 22;
    else if(atomic_num == 16) {
        type = 23;  // SA
        //if(atom.IsHbondDonor() || atom.MatchesSMARTS("[#16;!H0]"))  // SD
        if(atom.IsHbondDonor() || atom.ImplicitHydrogenCount()>0 || atom.ExplicitHydrogenCount()>0)
            type = 24;
        if(atom.MatchesSMARTS("[#16]=[#8]"))  // SO
            type = 25;
    }
    else if(atomic_num == 17)    // Cl
        type = 28;
    else if(atomic_num == 35)    // Br
        type = 29;
    else if(atomic_num == 53)    // I
        type = 30;
    else if(atomic_num == 30)    // Zn
        type = 31;
    else if(atomic_num == 25)    // Mn
        type = 31;
    else if(atomic_num == 12)    // Mg
        type = 31;
    else if(atomic_num == 26)    // Fe
        type = 31;
    else if(atomic_num == 23)    // V
        type = 31;
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
        FOR_NBORS_OF_ATOM(nbor,atom) {
            if(!nbor->IsCarbon() && !nbor->IsHydrogen())
                ++num_connected_het;
            if(isChargedOxygen(*nbor)) {   // CO
                type = 5;
                break;
            }
            if(isChargedNitrogen(*nbor)) { // CN
                type = 6;
                break;
            }
        }
        if(type == 0) {
            if(atom.IsAromatic()) {
                //if(atom.MatchesSMARTS("c~[!#6;!#1]"))   // cP
                if(num_connected_het)
                    type = 4;
                else                                    // cF
                    type = 3; 
            } else {
                //if(atom.MatchesSMARTS("[#6]~[!#6;!#1]"))  // CP
                if(num_connected_het)
                    type = 2;
                else                                      // CF
                    type = 1;
            }
        }
    }
    else if(atomic_num == 7) {
        type = 7;   // NA
        //if(atom.IsHbondAcceptor())  // NA
        //    type = 7;
        //if(atom.IsHbondDonor() || atom.MatchesSMARTS("[#7;!H0]"))  // ND
        if(atom.IsHbondDonor() || atom.ImplicitHydrogenCount()>0 || atom.ExplicitHydrogenCount()>0)
            type = 9;
        if(isChargedNitrogen(atom))  // NC
            type = 8;
    }
    else if(atomic_num == 8) {
        if(atom.MatchesSMARTS("[#8;H2]"))  // OW
            type = 13;
        else {
            type = 11;  // OA
            //if(atom.IsHbondAcceptor())     // OA
            //    type = 11;
            //if(atom.IsHbondDonor() || atom.MatchesSMARTS("[#8;!H0]")) // OD
            if(atom.IsHbondDonor() || atom.ImplicitHydrogenCount()>0 || atom.ExplicitHydrogenCount()>0)
                type = 12;
            if(isChargedOxygen(atom))      // OC
                type = 10;
        }
    }
    else if(atomic_num == 16) {
        type = 14;    // SA
        //if(atom.IsHbondAcceptor())   // SA
        //    type = 14;
        //if(atom.IsHbondDonor() || atom.MatchesSMARTS("[#16;!H0]"))  // SD
        if(atom.IsHbondDonor() || atom.ImplicitHydrogenCount()>0 || atom.ExplicitHydrogenCount()>0)
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

