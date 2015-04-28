/*=============================================================================
#     FileName: pmf_atom_pairs.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-08-15 09:14:44
#   LastChange: 2014-08-18 05:38:27
#      History:
=============================================================================*/
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include "atom_type.h"

using namespace std;
using namespace OpenBabel;

void get_protein_and_out_name(string &lig_name, string &pro_name, string &out_name);
void do_for_each(string &lig_name, string &pro_name, string &out_name);
double cutoff = 12.0;

#ifdef SEP
#undef SEP
#endif

#ifndef SEP
#ifdef WIN32
#define SEP "\\"
#else
#define SEP "/"
#endif
#endif

int main(int argc, char *argv[])
{
    if(argc!=3 && argc!=5 && argc!=7) {
        cerr << endl << "Usage: " << argv[0] << " [options]" << endl
			<< endl << "[options]" << endl
			<< "  -f in.pdb" << endl
			<< "  -l in.list" << endl
			<< "  --list list" << endl
            << endl << "Attention:" << endl
			<< "  1. if '-f' or '-l' is given, name of ligand should be like cmet_xxx_{pdbid}.pdb" << endl
            << "     and proteins will be read from /home/xmluo/jlpeng/cMet/pdb_protein/xxxx_protein.pdb." << endl
			<< "     atom pairs will be saved in file named 'cmet_xxx_{pdbid}.num'" << endl
			<< "  2. if '--list' is given, each line should be 'protein_name ligand_name'" << endl
			<< "     atom pairs will be saved in file named 'BASE{ligand_name}.num'" << endl
			<< "  3. *.num will be saved in the current directory" << endl
            << endl;
        exit(EXIT_FAILURE);
    }
    string mol_file("");
    string list_file("");
	string list("");
    for(int i=1; i<argc; ++i) {
        if(strcmp(argv[i], "-f") == 0)
            mol_file = argv[++i];
        else if(strcmp(argv[i], "-l") == 0)
            list_file = argv[++i];
		else if(strcmp(argv[i], "--list") == 0)
			list = argv[++i];
        else {
            cerr << "Error: invalid option " << argv[i] << endl;
            exit(EXIT_FAILURE);
        }
    }

	string lig_name(""), pro_name(""), out_name("");
    if(mol_file.size()) {
		get_protein_and_out_name(mol_file, pro_name, out_name);
        do_for_each(mol_file, pro_name, out_name);
	}
    if(list_file.size()) {
        ifstream inf(list_file.c_str());
        string line;
        while(getline(inf,line)) {
            if(line.size()==0 || line[0]=='#')
                continue;
			get_protein_and_out_name(line, pro_name, out_name);
            do_for_each(line, pro_name, out_name);
        }
        inf.close();
    }
	if(list.size()) {
		ifstream inf(list.c_str());
		string line;
		while(getline(inf, line)) {
			if(line.size()==0 || line[0]=='#')
				continue;
			istringstream is(line);
			is >> pro_name >> lig_name;
			string::size_type j = lig_name.rfind(SEP);
			if(j == string::npos)
				j = 0;
			else
				++j;
			string::size_type i = lig_name.rfind(".");
			out_name = lig_name.substr(j,i-j) + ".num";
			do_for_each(lig_name, pro_name, out_name);
		}
		inf.close();
	}

    return 0;
}

// lig_name: cmet_xxx_pdbid.pdb
// pro_name: pdbid_protein.pdb
// out_name: ${lig_name%format}num
void get_protein_and_out_name(string &lig_name, string &pro_name, string &out_name)
{
	string::size_type j = lig_name.rfind(SEP);
	if(j == string::npos)
		j = 0;
	else
		++j;
	string::size_type i = lig_name.rfind(".");
	string _format = lig_name.substr(i+1);
	string lig_basename = lig_name.substr(j,i-j);
	out_name = lig_basename+".num";
	i = lig_basename.rfind("_");
	string pdbid = lig_basename.substr(i+1);
	pro_name = "/home/xmluo/jlpeng/cMet/pdb_protein/"+pdbid+"_protein.pdb";
}

void do_for_each(string &lig_name, string &pro_name, string &out_name)
{
    string::size_type i = lig_name.rfind(".");
    string lig_format = lig_name.substr(i+1);
	i = pro_name.rfind(".");
	string pro_format = pro_name.substr(i+1);

    OBConversion conv;
    OBMol lig;
    OBMol pro;
    if(!conv.SetInFormat(lig_format.c_str()) || !conv.ReadFile(&lig, lig_name)) {
        cerr << "Error: failed to read ligand from " << lig_name << endl;
        return ;
    }
	if(!conv.SetInFormat(pro_format.c_str()) || !conv.ReadFile(&pro, pro_name)) {
		cerr << "Error: failed to read protein from " << pro_name << endl;
		return ;
	}

    ofstream outf(out_name.c_str());
    if(!outf) {
        cerr << "Error: failed to open " << out_name << endl;
        return ;
    }

    FOR_ATOMS_OF_MOL(latom, lig) {
        if(latom->IsHydrogen())
            continue;
        //cout << "for ligand " << latom->GetIdx() << endl;
        int ltype_idx = parse_ligand_atom_v1(*latom);
        if(ltype_idx == 0)
            continue;
        FOR_ATOMS_OF_MOL(patom, pro) {
            if(patom->IsHydrogen())
                continue;
            double dist = patom->GetDistance(&*latom);
            if(dist >= cutoff)
                continue;
            int ptype_idx = parse_protein_atom_v1(*patom);
            if(ptype_idx==0 || ptype_idx==13)
                continue;
            outf << LTYPE[ltype_idx-1] << "," << latom->GetIdx() << "," << PTYPE[ptype_idx-1]
                << "," << patom->GetIdx() << "," << dist << endl;
        }
    }

    outf.close();
}
