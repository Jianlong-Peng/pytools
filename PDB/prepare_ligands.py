'''
#=============================================================================
#     FileName: prepare_ligands.py
#         Desc: 
#       Author: 
#        Email: 
#     HomePage: 
#      Version: 0.0.1
#   LastChange: 2014-06-26 06:53:27
#      History:
#=============================================================================
'''
import os
import sys
import pybel


def main(argv=sys.argv):
	if len(argv) != 3:
		print ""
		print "  Usage: %s in_dir out_dir"%argv[0]
		print ""
		sys.exit(1)
	
	names = filter(lambda x: ("_lig" in x and x.endswith(".pdb")), os.listdir(argv[1]))
	for name in names:
		print "  >",name
		sys.stdout.flush()
		basename = os.path.basename(name)
		_format = basename[basename.rfind(".")+1:]
		basename = basename[:basename.rfind(".")]
		out_name = os.path.join(argv[2], basename+".sdf")
		if os.path.exists(out_name):
			continue
		mol = pybel.readfile(_format, name).next()
		mol.OBMol.CorrectForPH()
		mol.addh()
		mol.title = basename
		mol.write("sdf", out_name, overwrite=True)

main()

