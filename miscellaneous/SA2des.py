'''
#=============================================================================
#     FileName: SA2des.py
#         Desc: to create a descriptor from SAs
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-08-13 10:21:37
#   LastChange: 2012-08-13 10:42:08
#      History:
#=============================================================================
'''

import sys

def main(argv=sys.argv):
    if len(argv)!=4 and len(argv)!=5:
        print "Usage: SA2des.py want_frags mol_vs_frag result.txt [if_binary]"
        print "want_frags: SA name per line"
        print "mol_vs_frag: SAs matching information"
        print "result.txt: each line will be like \"mol value\""
        print "[if_binary]: <default: 1>"
        print " 1 - feature generated has binary value."
        print "     if mol matches at least one of SAs, it has value 1, otherwise 0."
        print " 2 - number of matched SAs as feature value."
        sys.exit(1)

    if len(argv) == 4:
        binary = 1
    else:
        binary = int(argv[4])
        if binary not in (1,2):
            print "[if_binary] should be either 1 or 2, but",binary,"given!"
            sys.exit(1)

    inf_want = open(argv[1],"r")
    wants = [line.strip() for line in inf_want]
    inf_want.close()

    inf_match = open(argv[2],"r")

    line = inf_match.readline()     #SA names
    line = line.split()
    try:
        indices = [line.index(item) for item in wants]
    except ValueError:
        print >>sys.err, "Sorry, but",argv[2],"doesn't have SA",item
        inf_match.close()
        sys.exit(1)
    
    outf = open(argv[3],"w")
    print >>outf, "mol SA_value"
    for line in inf_match:
        line = line.split()
        print >>outf, line[0],
        count = 0
        for i in indices:
            count += int(line[i])
        if binary == 1:
            print >>outf, (1 if count>=1 else 0)
        else:
            print >>outf, count

    inf_match.close()
    outf.close()

main()




