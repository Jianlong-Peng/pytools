/*=============================================================================
#     FileName: pmf.h
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-03-22 18:12:57
#   LastChange: 2014-06-19 15:38:24
#      History:
=============================================================================*/
#ifndef  PMF_H
#define  PMF_H

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <openbabel/mol.h>

const int Rmax = 12;
const double dr = 0.2;

std::map<std::string, std::vector<int> > count_number_of_pairs(
        OpenBabel::OBMol &protein, OpenBabel::OBMol &ligand,
        int *ptypes=NULL);

void write_number_of_pairs(
        std::map<std::string, std::vector<int> > &pmf_number,
        std::ostream &os);

#endif   /* ----- #ifndef PMF_H  ----- */

