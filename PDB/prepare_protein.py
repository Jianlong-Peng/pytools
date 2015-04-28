#Usage: ~/schrodinger/run prepare_protein.py in_dir out_dir

import sys
import os
from schrodinger import structure,structureutil
from schrodinger.application.prime import primefix


def main(argv=sys.argv):
    if len(argv)!=3 and len(argv)!=4:
        print "Usage: ~/schrodinger/run %s in_dir out_dir [finish_list]"%argv[0]
        sys.exit(1)
    
    finish_list = set([])
    if len(argv) == 4:
        inf = open(argv[3],"r")
        finish_list = [line.strip() for line in inf]
        inf.close()
    names = filter(lambda x: "_pro" in x and x.endswith(".pdb"), os.listdir(argv[1]))
    for name in names:
        if name in finish_list:
            continue
        print "  >",name
        sys.stdout.flush()
        out_name = os.path.join(argv[2],name)
        #basename = name[:name.rfind(".")]
        #out_name = os.path.join(argv[2],basename+"_prep.pdb")
        #read protein
        pro = structure.StructureReader(os.path.join(argv[1],name)).next()
        #fix using prime
        try:
            primefix.fix_for_prime(pro)
        except:
            print "primefix failed"
            print sys.exc_info()
        #re-add hydrogens
        structureutil.delete_hydrogens(pro)
        structureutil.add_hydrogens(pro)
        #save
        s = structure.Structure(pro.handle)
        s.write(out_name)
        #to check if protein has missing atom(s)
        """
        resnum = None
        for i in xrange(len(pro.atom)):
            if "i_m_pdb_convert_problem" in pro.atom[i+1].property:
                resnum = pro.atom[i+1].resnum
                break
        if resnum is not None:
            print "   ",name,"has i_m_pdb_covert_problem for residue",resnum
        """

main()

