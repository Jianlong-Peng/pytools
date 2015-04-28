/*=============================================================================
#     FileName: pmf.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-03-22 18:19:16
#   LastChange: 2014-06-20 10:04:30
#      History:
=============================================================================*/
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <utility>
#include <openbabel/mol.h>
#include "pmf.h"
#include "atom_type.h"

using namespace std;
using namespace OpenBabel;

map<string, vector<int> > count_number_of_pairs(
        OBMol &protein, OBMol &ligand, int *ptypes)
{

    int *ltypes = parse_ligand_v1(ligand);

    map<string, vector<int> > pmf_number;
    int num_bins = static_cast<int>(Rmax / dr);
    for(vector<string>::iterator i=PTYPE.begin(); i!=PTYPE.end(); ++i)
        for(vector<string>::iterator j=LTYPE.begin(); j!=LTYPE.end(); ++j)
            pmf_number.insert(make_pair(*i+"_"+*j, vector<int>(num_bins, 0)));

    FOR_ATOMS_OF_MOL(patom, protein) {
        unsigned pidx = patom->GetIdx() - 1;
        int ptype_idx;
        if(ptypes == NULL)
            ptype_idx = parse_protein_atom_v1(*patom);
        else
            ptype_idx = ptypes[pidx];
        if(ptype_idx == 0)
            continue;
        FOR_ATOMS_OF_MOL(latom, ligand) {
            unsigned lidx = latom->GetIdx() - 1;
            int ltype_idx = ltypes[lidx];
            if(ltype_idx == 0)
                continue;
            double dist = patom->GetDistance(&*latom);
            if(dist > Rmax)
                continue;
            string key = PTYPE[ptype_idx-1] + "_" + LTYPE[ltype_idx-1];
            int index = static_cast<int>(ceil(dist / dr)) - 1;
            pmf_number[key][index] += 1;  // key must be in pmf_number !!!
        }
    }

    free(ltypes);

    return pmf_number;
}

void write_number_of_pairs(map<string, vector<int> > &pmf_number, ostream &os)
{
    for(map<string,vector<int> >::iterator i=pmf_number.begin(); i!=pmf_number.end(); ++i) {
        os << i->first;
        for(vector<int>::iterator j=(i->second).begin(); j!=(i->second).end(); ++j)
            os << " " << *j;
        os << endl;
    }
}

