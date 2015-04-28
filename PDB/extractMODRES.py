'''
#=============================================================================
#     FileName: extractMODRES.py
#         Desc: 
#       Author: 
#        Email: 
#     HomePage: 
#      Version: 0.0.1
#   LastChange: 2014-06-05 10:36:20
#      History:
#=============================================================================
'''
import sys

def main(argv=sys.argv):
    if len(argv) != 2:
        print "\n  Usage: %s in_list output"%argv[0]
        print "  in_list: pdb file in each line"
        print "  output : where modified resdiues will be saved"
        print ""
        sys.exit(1)
    
    inf = open(argv[1],"r")
    for line in inf:
        inf_each = open(line.strip(),"r")
        for line_each in inf_each:
            if line_each.startswith("MODRES"):
               print line_each[12:15],line_each[24:27],line.strip()
        inf_each.close()
    inf.close()
    """
    residues = {}
    inf = open(argv[1],"r")
    for line in inf:
        inf_each = open(line.strip(),"r")
        for line_each in inf_each:
            if line_each.startswith("MODRES"):
               new_name = line_each[12:15]
               old_name = line[24:27]
               name = new_name+"_"+old_name
               if not residues.has_key(name):
                   residues[name] = 0
               residues[name] += 1
        inf_each.close()
    inf.close()
    
    outf = open(argv[2],"w")
    print >>outf, "new old frequency"
    for key,value in residues.iteritems():
        new_name,old_name = key.split("_")
        print >>outf, new_name,old_name,value
    outf.close()
    """

main()

