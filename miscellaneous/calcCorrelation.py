'''
#=============================================================================
#     FileName: calcCorrelation.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-12-05 10:21:22
#   LastChange: 2013-12-05 10:31:54
#      History:
#=============================================================================
'''
import sys
import numpy as np

def main(argv=sys.argv):
    if len(argv)!=3 and len(argv)!=4:
        print "\n  Usage: %s input output [sep]"%argv[0]
        print "  [sep]: <default: space>\n"
        sys.exit(1)

    if len(argv) == 4:
        sep = argv[3]
    else:
        sep = ""

    X = []
    y = []
    inf = open(argv[1],"r")
    line = inf.readline()
    if sep == "":
        names = line.split()[1:]
    else:
        names = line.strip().split(sep)[1:]
    for line in inf:
        if sep == "":
            line = line.split()
        else:
            line = line.strip().split(sep)
        y.append(float(line[0]))
        X.append(map(float,line[1:]))
    inf.close()

    X = np.asarray(X).T
    y = np.asarray(y)

    outf = open(argv[2],"w")
    print >>outf, "desc\tr"
    for i in xrange(X.shape[0]):
        r = np.corrcoef(X[i],y)
        print >>outf, "%s\t%g"%(names[i],r[0][1])
    outf.close()

main()

