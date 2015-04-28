'''
#=============================================================================
#     FileName: calc_potential.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-07-22 10:30:24
#   LastChange: 2014-07-22 10:32:54
#      History:
#=============================================================================
'''
import sys
import pickle
import dist


def main(argv=sys.argv):
    if len(argv) != 3:
        print "\n  Usage: %s num_pairs.dat potential.dat"%argv[0]
        print ""
        sys.exit(1)

    inf = open(argv[1],"r")
    dc = pickle.load(inf)
    inf.close()

    p = dist.PTS(dc)
    outf = open(argv[2],"w")
    pickle.dump(p, outf)
    outf.close()

main()

