'''
#=============================================================================
#     FileName: pmf_atom_pairs.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-08-14 14:22:43
#   LastChange: 2014-08-15 01:22:37
#      History:
#=============================================================================
'''
import sys
from glob import glob
import pybel
import pmf_atom_typer

cutoff = 12.0

pro_atom = {'3Q6W':[], '3U6H':[], '3VW8':[], '3ZCL':[]}

def calc_pro_types(pdb_id):
    name = "/home/xmluo/jlpeng/cMet/pdb_protein/%s_protein.pdb"%pdb_id
    pro = pybel.readfile("pdb", name).next()
    pro_atom[pdb_id] = [None]*len(pro.atoms)
    for i in xrange(len(pro.atoms)):
        pro_atom[pdb_id][i] = pmf_atom_typer.ProAtomTyper(pro.atoms[i].OBAtom)

def do_for_each(name):
    _format = name[name.rfind(".")+1:]
    basename = name[:name.rfind(".")]
    out_name = basename + ".num"
    pdb_id = basename.split("_")[-1]
    #if len(pro_atom[pdb_id]) == 0:
    #    print "to do atom typing for %s..."%pdb_id,
    #    sys.stdout.flush()
    #    calc_pro_types(pdb_id)
    #    print "  Done"
    #    sys.stdout.flush()
    pdb_name = "/home/xmluo/jlpeng/cMet/pdb_protein/%s_protein.pdb"%pdb_id
    lig = pybel.readfile(_format, name).next()
    pro = pybel.readfile("pdb", pdb_name).next()
    outf = open(out_name, "w")
    for al in lig.atoms:
        if al.OBAtom.IsHydrogen():
            continue
        ltype = pmf_atom_typer.LigAtomTyper(al.OBAtom)
        lidx = al.OBAtom.GetIdx()
        if ltype == "":
            continue
        for j in xrange(len(pro.atoms)):
            ap = pro.atoms[j]
            if ap.OBAtom.IsHydrogen():
                continue
            dist = ap.OBAtom.GetDistance(al.OBAtom)
            if dist >= cutoff:
                continue
            #ptype = pro_atom[pdb_id][j]
            ptype = pmf_atom_typer.ProAtomTyper(ap.OBAtom)
            pidx = ap.OBAtom.GetIdx()
            if ptype in ('','OW'):
                continue
            print >>outf, "%s,%d,%s,%d,%.6f"%(ltype, lidx, ptype, pidx, dist)
    outf.close()


def main(argv=sys.argv):
    if len(argv) < 2:
        print "\n  Usage: %s mol_in[...]"%argv[0]
        print "  mol_in: should be like cmet_xxx_{pdbid}.pdb"
        print "\n  information of atom pair for each molecule"
        print "  will be saved in a file named basename{mol_in}.num"
        print ""
        sys.exit(1)

    names = []
    for name in argv[1:]:
        for item in glob(name):
            names.append(item)
    print names
    for name in names:
        do_for_each(name)

main()

