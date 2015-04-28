'''
#=============================================================================
#     FileName: csv2arff.py
#         Desc: convert a csv format into .arff format
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-09-15 14:11:01
#   LastChange: 2012-09-15 14:37:29
#      History:
#=============================================================================
'''

import sys
import os

def main(argv=sys.argv):
    if len(argv)!=3 and len(argv)!=4:
        print "Usage: csv2arff.py in.csv out.arff [feature_def]"
        print "- in.csv: 1st row with feature names, and 1st col with class labels or ys"
        print "- outf.arff"
        print "- [feature_def] - each line will be: feature_name continue[,discrete]"
        print "  if not given, all attributes will be cheated as continue"
        sys.exit(1)

    inf = open(argv[1],"r")
    outf = open(argv[2],"w")
    print >>outf, "@relation",os.path.basename(argv[1])
    line = inf.readline()
    feature_names = line[:-1].split(",")[1:]
    feature_types = ["real" for _ in xrange(len(feature_names))]

    if len(argv) == 4:
        inf_des = open(argv[3],"r")
        for _line in inf_des:
            _line = _line.split()
            i = feature_names.index(_line[0])
            feature_types[i] = _line[1]
        inf_des.close()
    
    for i in xrange(len(feature_names)):
        print >>outf, "@attribute",feature_names[i],feature_types[i]

    ys = []
    dataset = []
    for line in inf:
        i = line.index(",")
        ys.append(line[:i])
        dataset.append(line[i+1:-1])
    tmp_ys = set(ys)
    if len(tmp_ys) < len(ys)/3.:   #naive??
        print >>outf, "@attribute class {"+(",".join(tmp_ys))+"}"
    else:
        print >>outf, "@attribute class real"
    print >>outf, "\n@DATA"

    for i in xrange(len(dataset)):
        print >>outf, dataset[i]+","+ys[i]

    inf.close()
    outf.close()

main()


