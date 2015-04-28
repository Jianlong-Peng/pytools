'''
#=============================================================================
#     FileName: svm2csv.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-12-18 16:46:01
#   LastChange: 2012-12-18 17:04:06
#      History:
#=============================================================================
'''

import sys

def main(argv=sys.argv):
    if len(argv) != 3:
        print "Usage: svm2csv.py in.svm out.csv"
        sys.exit(1)

    tags = []
    inf = open(argv[1],"r")
    for line in inf:
        line = line.split()
        tags.extend([int(item.split(":")[0]) for item in line[1:]])
    max_tag = max(tags)

    inf.seek(0)
    outf = open(argv[2],"w")
    outf.write("label")
    for i in xrange(max_tag):
        outf.write(" %d"%(i+1))
    outf.write("\n")
    for line in inf:
        line = line.split()
        outf.write(line[0])
        i = 0
        for item in line[1:]:
            i += 1
            index,value = item.split(":")
            while i < int(index):
                outf.write(" 0")
                i += 1
            outf.write(" "+value)
        while i < max_tag:
            outf.write(" 0")
            i += 1
        outf.write("\n")
        

    inf.close()
    outf.close()

main()


