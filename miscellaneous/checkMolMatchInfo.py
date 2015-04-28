'''
#=============================================================================
#     FileName: checkMolMatchInfo.py
#         Desc: to get info about wether mol matches any SAs
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-08-16 13:47:14
#   LastChange: 2012-08-16 14:01:13
#      History:
#=============================================================================
'''

import sys
import os

def main(argv=sys.argv):
    if len(argv) != 4:
        print "Usage: checkMolMatchInfo.py mol_vs_frag want_frag result.txt"
        sys.exit(1)

    assert os.path.isfile(argv[1]) and os.path.isfile(argv[2])

    inf_want_frag = open(argv[2],"r")
    want_frags = [line.strip() for line in inf_want_frag]
    inf_want_frag.close()

    inf = open(argv[1],"r")
    line = inf.readline()
    line = line.split()[1:]
    sa_id = [i for i in xrange(len(line)) if want_frags.count(line[i])]
    outf = open(argv[3],"w")
    print >>outf, "mol numOfSAs_contains"
    for line in inf:
        line = line.split()
        print >>outf, line[0],
        line = [int(item) for item in line[1:]]
        count = 0
        for _id in sa_id:
            count += line[_id]
        print >>outf, count

    inf.close()
    outf.close()

main()


