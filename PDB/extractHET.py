'''
#=============================================================================
#     FileName: extractHET.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-06-04 19:41:50
#   LastChange: 2014-06-04 11:55:36
#      History:
#=============================================================================
'''
import sys


def main(argv=sys.argv):
    if len(argv) != 2:
        print "\n  Usage: %s in_list output"%argv[0]
	print ""
        print "  in_list: pdb file in each line"
        print "  output : where HET names will be saved"
        print ""
        sys.exit(1)

    inf = open(argv[1],"r")
    for name in inf:
        inf_each = open(name.strip(), "r")
        for line in inf_each:
            if line.startswith("HET   "):
                print line[7:10]+"\t"+name.strip()
        inf_each.close()
    inf.close()

    """
    het_names = {}

    inf = open(argv[1],"r")
    for name in inf:
        inf_each = open(name.strip(),"r")
        for line in inf_each:
            if line.startswith("HET   "):
                het = line[7:10]
                if not het_names.has_key(het):
                    het_names[het] = 0
                het_names[het] += 1
        inf_each.close()
    inf.close()
    print "totally found %d unique HET names"%len(het_names.keys())

    outf = open(argv[2],"w")
    print >>outf, "HET frequency"
    for key,value in het_names.iteritems():
        print >>outf, key,value
    outf.close()
    """

main()


