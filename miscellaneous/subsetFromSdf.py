'''
#=============================================================================
#     FileName: subsetFromSdf.py
#         Desc: to get a subset of molecules from multi-molecule sdf file
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-09-27 09:57:44
#   LastChange: 2012-09-27 10:06:34
#      History:
#=============================================================================
'''

import sys

def main(argv=sys.argv):
    if len(argv) != 4:
        print "Usage: subsetFromSdf.py in.sdf want_names.txt out.sdf"
        print "want_names.txt: a name per line"
        print "recommend: `in.sdf` should not be too large!!"
        sys.exit(1)

    inf = open(argv[1],"r")
    mols = {}
    line = inf.readline()
    content = ""
    while line != "":
        if line.strip() == "":
            line = inf.readline()
            continue
        name = line.strip()
        content = line
        line = inf.readline()
        while line.strip() != "$$$$":
            content += line
            line = inf.readline()
        content += line.strip()
        mols[name] = content
        line = inf.readline()
    inf.close()

    inf_want = open(argv[2],"r")
    outf = open(argv[3],"w")
    for name in inf_want:
        name = name.strip()
        print >>outf, mols[name]

    inf_want.close()
    outf.close()

main()



