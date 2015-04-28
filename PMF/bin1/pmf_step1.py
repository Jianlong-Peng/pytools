'''
#=============================================================================
#     FileName: pmf_step1.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-03-22 15:51:01
#   LastChange: 2014-07-22 15:25:25
#      History:
#=============================================================================
'''
import sys
import os
from getopt import getopt
import pybel
import pmf
import atom_type

def exit_with_help(name):
    print "\nOBJ:"
    print "  to count numbers of atom pairs ij for each complex, n_{ij,m}(r)"
    print "\nUsage:"
    print "  %s [options] outdir"%name
    print "\n[options]"
    print "  -p protein: protein file, (pdb)"
    print "  -l ligand: ligand file, (pdb, mol, sdf, mol2 etc.)"
    print "  -c complex: pdb file containg protein and ligand"
    print "  --ligand-list file: file contains ligand files"
    print "           in this case, '-p protein' must be given!!"
    print "  --list file: each line should be 'complex' or 'protein ligand'"
    print "  --type n: to specify the atom type being used <default: 0>"
    print "    0 - PMF04 atom type"
    print "  --show-types: to display atom types to be used"
    print "\nAttention:"
    print "  1. the numbers of atom pairs i of each complex will be saved separately."
    print "     the file wil be named as {basename}_num.txt"
    print "  2. it's useful to use '--ligand-list' when number of atom pairs is to be"
    print "     calculated for many ligands against one protein"
    print ""
    sys.exit(1)

def main(argv=sys.argv):
    if len(argv) == 1:
        exit_with_help(argv[0])

    options,args = getopt(argv[1:], "p:l:c:",["list=","type=","show-types","ligand-list="])
    protein_file = None
    ligand_file = None
    ligand_list = None
    complex_file = None
    list_file = None
    _type = 0
    show = False
    for opt,value in options:
        if opt == "-p":
            protein_file = value
        elif opt == "-l":
            ligand_file = value
        elif opt == "-c":
            complex_file = value
        elif opt == "--list":
            list_file = value
        elif opt == "--type":
            _type = int(value)
            assert _type in atom_type.available_types
        elif opt == "--show-types":
            show = True
        elif opt == "--ligand-list":
            ligand_list = val
        else:
            print >>sys.stderr, "Error: invalid option",opt
            sys.exit(1)

    if show:
        print atom_type.atom_type_description[_type]
        return 0

    if (ligand_list is not None) and (protein_file is None):
        print >>sys.stderr, "Error: '--ligand-list' is given without '-p'"
    if (protein_file is not None) and (ligand_file is None and ligand_list is None):
        print >>sys.stderr, "Error: -p is given without -l or --ligand-list"
        sys.exit(1)
    if (protein_file is None) and (ligand_file is not None):
        print >>sys.stderr, "Error: -l is given without -p"
        sys.exit(1)

    assert len(args) == 1
    out_dir = args[0]
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    if ligand_file is not None:
        print "to process",protein_file,"&",ligand_file
        out_name = get_out_name(protein_file, ligand_file)
        out_name = os.path.join(out_dir, out_name+"_num.txt")
        assert not os.path.exists(out_name)
        outf = open(out_name, "w")
        print >>outf, "atom type %d"%_type
        #print >>outf, "Rmax=%g, dr=%g"%(pmf.Rmax, pmf.dr)
        count_protein_ligand(protein_file, ligand_file, outf, _type)
        outf.close()

    if ligand_list is not None:
        inf = open(ligand_list,"r")
        for line in inf:
            line = line.strip()
            if line=="" or line.startsiwht("#"):
                continue
            print "to process",protein_file,"&",line
            out_name = get_out_name(protein_file, line)
            out_name = os.path.join(out_dir, out_name+"_num.txt")
            assert not os.path.exists(out_name)
            outf = open(out_name,"w")
            print >>outf, "atom type %d"%_type
            count_protein_ligand(protein_file, line, outf, _type)
            outf.close()
        inf.close()

    if complex_file is not None:
        print >>sys.stderr, "Warning: '-c' is currently not supported!"
    
    if list_file is not None:
        inf = open(list_file,"r")
        for line in inf:
            if line.strip()=="" or line.startswith("#"):
                continue
            line = line.split()
            if len(line) == 1:
                print >>sys.stderr, "Warning: currently pdb file with complex is not supported!"
                print line.strip()
                continue
            elif len(line) == 2:
                print "to process",line[0],"&",line[1]
                out_name = get_out_name(line[0], line[1])
                out_name = os.path.join(out_dir, out_name+"_num.txt")
                assert not os.path.exists(out_name)
                outf = open(out_name,"w")
                print >>outf, "atom type %d"%_type
                #print >>outf, "Rmax=%g, dr=%g"%(pmf.Rmax, pmf.dr)
                count_protein_ligand(line[0],line[1],outf,_type)
                outf.close()
            else:
                print >>sys.stderr, "Error: each line of 'list_file' should be one or two file names"
                print >>sys.stderr, "but got '%s'"%(line.strip())
        inf.close()
    

def get_out_name(pro_file, lig_file):
    pro_basename = os.path.basename(pro_file)
    temp = pro_basename.split("_")
    if len(temp) == 1:
        pro_name = pro_basename[:pro_basename.rfind(".")]
    elif len(temp) == 2:   #pdb1cea_pro.pdb
        pro_name = temp[0]
    elif len(temp) == 3:  #pdb1ce8_pro_1.pdb
        pro_suffix = temp[2].split(".")[0]
        pro_name = temp[0]+"_"+pro_suffix
    else:
        raise ValueError(pro_file)

    lig_basename = os.path.basename(lig_file)
    temp = lig_basename.split("_")
    if len(temp) == 1:
        lig_name = lig_basename[:lig_basename.rfind(".")]
    elif len(temp) == 2:
        lig_name = temp[0]
    elif len(temp) == 3:
        lig_suffix = temp[2].split(".")[0]
        lig_name = temp[0]+"_"+lig_suffix
    else:
        lig_name = temp[0]+"_"+temp[1]
        #raise ValueError(lig_file)

    if pro_name == lig_name:
        return pro_name
    else:
        return pro_name+"_"+lig_name


def count_protein_ligand(protein_file, ligand_file, outf, _type):
    #read protein
    pro_format = protein_file[protein_file.rfind(".")+1:]
    pro_base_name = os.path.basename(protein_file)
    protein = pybel.readfile(pro_format, protein_file).next()
    ptypes = atom_type.parse_protein(protein.OBMol, _type)
    #read ligand and count number of atom pairs
    lig_format = ligand_file[ligand_file.rfind(".")+1:]
    lig_base_name = os.path.basename(ligand_file)
    inf_ligand = pybel.readfile(lig_format, ligand_file)
    i = 0
    for ligand in inf_ligand:
        i += 1
        print "  to count numbers of atom pairs for %dth ligand"%i
        ltypes = atom_type.parse_ligand(ligand.OBMol, _type)
        ap_pl = pmf.gen_atom_pair_PL(protein.OBMol, ligand.OBMol, _type, ptypes, ltypes)
        ap_ll = pmf.gen_atom_pair_LL(ligand.OBMol, _type, ltypes)
        #pmf_number = pmf.count_number_of_pairs(protein.OBMol, ligand.OBMol, ptypes)
        print >>outf, "> protein=%s; ligand=%s (%dth)"%(pro_base_name,lig_base_name,i)
        outf.writelines(ap_pl)
        outf.writelines(ap_ll)
        #pmf.write_number_of_pairs(pmf_number, outf)


main()

