/*=============================================================================
#     FileName: pmf_step1.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-03-22 18:40:17
#   LastChange: 2014-06-20 10:41:07
#      History:
=============================================================================*/
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cctype>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include "atom_type.h"
#include "pmf.h"


using namespace std;
using namespace OpenBabel;

bool count_protein_ligand(string &protein_file, string &ligand_file, ofstream &outf, bool verbose);

vector<string> split(string &line)
{
    vector<string> result;
    unsigned i = 0;
    unsigned j;
    while(i < line.size()) {
        while(i<line.size() && isspace(line[i]))
            ++i;
        if(i == line.size())
            break;
        j = i;
        while(j<line.size() && isspace(line[j]))
            ++j;
        result.push_back(line.substr(i,j-i));
        i = j;
    }
    return result;
}

#if defined(WIN32)
#define SEP "\\"
#else
#define SEP "/"
#endif

#define TEST_PMF

#if defined(TEST_PMF)
void exit_with_help(const char *name)
{
    cerr << endl << "OBJ" << endl
        << "  to count numbers of atom pairs ij for each complex m, n_ijm(r)" << endl
        << endl << "Usage" << endl
        << "  " << name << " [options] outfile" << endl
        << endl << "[options]" << endl
        << "  -p protein: protein file, (pdb)" << endl
        << "  -l ligand : ligand file, (pdb, mol, sdf, mol2 etc.)" << endl
        << "  -c complex: pdb file containing protein and ligand" << endl
        << "  --list file: each line should be 'complex' or 'protein ligand'" << endl
        << "  --verbose: if given, to display which complex is being processed" << endl
        << endl << "Attention" << endl
        << "  1. reference sphere with a radius of R=" << Rmax << endl
        << "  2. bin size: " << dr << endl
        << endl;
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    if(argc<4 || argc>11)
        exit_with_help(argv[0]);

    string protein_file("");
    string ligand_file("");
    string complex_file("");
    string list_file("");
    bool verbose(false);
    int i=1;
    for(; i<argc; ++i) {
        if(argv[i][0] != '-')
            break;
        if(strcmp(argv[i],"-p") == 0)
            protein_file = argv[++i];
        else if(strcmp(argv[i],"-l") == 0)
            ligand_file = argv[++i];
        else if(strcmp(argv[i],"-c") == 0)
            complex_file = argv[++i];
        else if(strcmp(argv[i],"--list") == 0)
            list_file = argv[++i];
        else if(strcmp(argv[i],"--verbose") == 0)
            verbose = true;
        else {
            cerr << "Error: invalid option " << argv[i] << endl;
            exit(EXIT_FAILURE);
        }
    }
    if(argc-i != 1) {
        cerr << "Error: invalid number of arguments" << endl;
        exit(EXIT_FAILURE);
    }

    if(protein_file.empty() && !ligand_file.empty()) {
        cerr << "Error: -p is given without -l" << endl;
        exit(EXIT_FAILURE);
    }
    if(!protein_file.empty() && ligand_file.empty()) {
        cerr << "Error: -l is given without -p" << endl;
        exit(EXIT_FAILURE);
    }

    ofstream outf(argv[i]);
    if(!outf) {
        cerr << "Error: failed to open file " << argv[i] << endl;
        exit(EXIT_FAILURE);
    }
    outf << "Rmax=" << Rmax << ", dr=" << dr << endl;
    //int num_bins = static_cast<int>(Rmax / dr);
    /*
    outf << "atom_pairs(pro_lig):";
    for(vector<string>::iterator i=PTYPE.begin(); i!=PTYPE.end(); ++i)
        for(vector<string>::iterator j=LTYPE.begin(); j!=LTYPE.end(); ++j)
            outf << " " << *i << "_" << *j;
    outf << endl;
    */
    // 1. protein & ligand
    if(!protein_file.empty() && !ligand_file.empty()) {
        if(verbose)
            cout << "to process " << protein_file << " & " << ligand_file << "..." << endl;
        bool success = count_protein_ligand(protein_file, ligand_file, outf, verbose);
        if(verbose)
            cout << "=>" << (success?"success":"fail") << endl << endl;
    }

    // 2. complex
    if(!complex_file.empty())
        cerr << "Warning: '-c' is currently not supported" << endl;

    // 3. list
    if(!list_file.empty()) {
        ifstream inf(list_file.c_str());
        if(!inf) {
            cerr << "Error: failed to open file " << list_file << endl;
            outf.close();
            exit(EXIT_FAILURE);
        }
        string line;
        int line_no = 0;
        while(getline(inf,line)) {
            ++line_no;
            if(line.size()==0 || line[0]=='#')
                continue;
            vector<string> temp = split(line);
            if(temp.size() == 1) {
                cout << "Warning: only one file found in Line "<< line_no << endl
                    << "         pdb file with complex is currently not supported" << endl;
                continue;
            }
            else if(temp.size() == 2) {
                if(verbose)
                    cout << "to process " << temp[0] << " & " << temp[1] << "..." << endl;
                bool success = count_protein_ligand(temp[0], temp[1], outf, verbose);
                if(verbose)
                    cout << "=>" << (success?"success":"fail") << endl << endl;
            }
            else {
                cout << "Warning: more than 2 files were found in Line " << line_no << endl;
                continue;
            }
        }
        inf.close();
    }

    outf.close();

    return 0;
}
#else
int main(int argc, char *argv[])
{
    if(argc != 3) {
        cerr << endl << "OBJ" << endl
            << "  to count number of atom pairs ij for each complex m, n_ijm(r)" << endl
            << endl << "Usage" << endl
            << "  " << argv[0] << " in_list out_dir" << endl
            << endl << "Arguments" << endl
            << "  in_list: each line should be 'dir pdb_id' or 'pdb_id'" << endl
            << "  out_dir: where to save files containing number of atom pairs" << endl
            << endl;
        exit(EXIT_FAILURE);
    }

    ifstream inf(argv[1]);
    if(!inf) {
        cerr << "Error: failed to open file " << argv[1] << endl;
        exit(EXIT_FAILURE);
    }
    string line;
    string pro_file, lig_file, out_file;
    int line_no = 0;
    while(getline(inf,line)) {
        ++line_no;
        if(line.size()==0 || line[0]=='#')
            continue;
        vector<string> temp = split(line);
        if(temp.size() == 1) {
            pro_file = temp[0]+"_pro.pdb";
            lig_file = temp[0]+"_lig.pdb";
            out_file = string(argv[2])+SEP+temp[0]+"_num.txt";
        }
        else if(temp.size() == 2) {
            pro_file = temp[0]+SEP+temp[1]+"_pro.pdb";
            lig_file = temp[0]+SEP+temp[1]+"_lig.pdb";
            out_file = string(argv[2])+SEP+temp[1]+"_num.txt";
        }
        else {
            cout << "Warning: line " << line_no << " contains more than 2 items" << endl;
            continue;
        }
        cout << "to process " << pro_file << " & " << lig_file << endl;
        ofstream outf(out_file);
        if(!outf) {
            cout << "  Error: failed to open output file " << out_file;
            continue;
        }
        outf << "Rmax=" << Rmax << ", dr=" << dr << endl;
        //int num_bins = static_cast<int>(Rmax / dr);
        /*
        outf << "atom_pairs(pro_lig):";
        for(vector<string>::iterator i=PTYPE.begin(); i!=PTYPE.end(); ++i)
            for(vector<string>::iterator j=LTYPE.begin(); j!=LTYPE.end(); ++j)
                outf << " " << *i << "_" << *j;
        outf << endl;
        */
        bool success = count_protein_ligand(pro_file, lig_file, outf, true);
        cout << "=>success?" << (success?"true":"false") << endl;
        outf.close();
    }
    inf.close();

    return 0;
}
#endif

bool count_protein_ligand(string &protein_file, string &ligand_file, ofstream &outf, bool verbose)
{
    // read protein
    string::size_type i = protein_file.rfind(".");
    if(i == string::npos) {
        cout << "  Error: can't find the format of " << protein_file << endl;
        return false;
    }
    string pro_format = protein_file.substr(i+1);
    OBConversion pro_conv;
    OBMol protein;
    if(!(pro_conv.SetInFormat(pro_format.c_str()) && pro_conv.ReadFile(&protein, protein_file))) {
        cout << "  Error: failed to read protein from " << protein_file << endl;
        return false;
    }
    if(verbose)
        cout << "  to parse protein atom type" << endl;
    int *ptypes = parse_protein_v1(protein);
    // read ligand and count number of atom pairs
    i = ligand_file.rfind(".");
    if(i == string::npos) {
        cout << "  Error: can't find the format of " << ligand_file << endl;
        return false;
    }
    string lig_format = ligand_file.substr(i+1);
    OBConversion lig_conv;
    OBMol ligand;
    if(!lig_conv.SetInFormat(lig_format.c_str())) {
        cout << "  Error: failed to set input format " << lig_format << " of ligand file" << endl;
        return false;
    }
    bool notatend = lig_conv.ReadFile(&ligand, ligand_file);
    int j = 0;
    while(notatend) {
        ++j;
        if(verbose)
            cout << "  to count number of atom pairs for " << j << "th ligand" << endl;
        map<string,vector<int> > pmf_number = count_number_of_pairs(protein, ligand, ptypes);
        outf << "> protein=" << protein_file << "; ligand=" << ligand_file << "_" << j << endl;
        write_number_of_pairs(pmf_number, outf);
        notatend = lig_conv.Read(&ligand);
    }

    free(ptypes);

    return true;
}

