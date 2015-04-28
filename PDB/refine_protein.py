'''
#=============================================================================
#     FileName: refine_protein.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-06-28 14:03:58
#   LastChange: 2014-06-29 08:07:01
#      History:
#=============================================================================
'''
import sys
import os
from schrodinger import structure

inp_info = """JOB_TYPE\tREFINE
PRIME_TYPE\tSTDE_PRED
SELECT\tpick
NUM_SC_OUTPUT_STRUCT\t1
USE_CRYSTAL_SYMMETRY\tno
USE_RANDOM_SEED\tno
SEED\t0
EXT_DIEL\t80.00
USE_MEMBRANE\tno
HOST\tlocalhost
"""

def main(argv=sys.argv):
    if len(argv) != 2:
        print "\n  Usage: ~/schrodinger/run %s check_protein.log"%argv[0]
        print "  OBJ: parse 'check_protein.log', and generate .inp file for 'refinestruct'"
        print ""
        sys.exit(1)

    root = ""
    inf = open(argv[1],"r")
    line = inf.readline()
    while line != "":
        if line.startswith("to check"):
            root = line.split()[-1]
            line = inf.readline()
            continue
        if line.startswith("  >"):
            res_chain = set([])
            pdb_name = line.split()[-1]
            print "%s/%s"%(root,pdb_name)
            sys.stdout.flush()
            line = inf.readline()
            while line!="" and not line.startswith("to check") and not line.startswith("  >"):
                if "missing atom" in line:
                    temp = line.split(",")
                    res = temp[0].split("residue")[-1].strip()
                    #assert res.isdigit()
                    chain = temp[1].split("'")[-2]
                    res_chain.add(chain+":"+res)
                line = inf.readline()
            res_chain = list(res_chain)
            basename = pdb_name[:pdb_name.rfind(".")]
            if (len(res_chain)!=0) and (not os.path.exists(os.path.join(root, basename+".inp"))):
                #generate {base}.mae
                pro = structure.StructureReader(os.path.join(root,pdb_name)).next()
                s = structure.Structure(pro.handle)
                s.write(os.path.join(root, basename+".mae"))
                #generate {bsae}.inp
                outf = open(os.path.join(root, basename+".inp"), "w")
                print >>outf, "STRUCT_FILE\t%s"%(basename+".mae")
                print >>outf, "JOB_TYPE\tREFINE\nPRIME_TYPE\tSIDE_PRED\nSELECT\tpick"
                for i in xrange(len(res_chain)):
                    print >>outf, "RESIDUE_%d\t%s"%(i, res_chain[i])
                print >>outf, """NUM_SC_OUTPUT_STRUCT\t1
USE_CRYSTAL_SYMMETRY\tno
USE_RANDOM_SEED\tno
SEED\t0
EXT_DIEL\t80.00
USE_MEMBRANE\tno
HOST\tlocalhost
"""
                outf.close()
        else:
            line = inf.readline()

    inf.close()

main()

