'''
#=============================================================================
#     FileName: genSVMformat.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-12-05 14:44:22
#   LastChange: 2013-12-05 18:25:29
#      History:
#=============================================================================
'''
import sys
from getopt import getopt
import numpy as np

def main(argv=sys.argv):
    if len(argv) not in (3,5,7):
        print "\n    OBJ: to generate LIBSVM format"
        print "         suppose that y-values were given in 1st column"
        print "\n  Usage: %s [options] input output"%argv[0]
        print "  [options]"
        print "  --sep delimiter   <default: white space>"
        print "  --mark file"
        print "         if given, a file contains information about"
        print "         descriptors to be used to construct the model."
        print "         each line should be a descriptor name to be used!"
        print "         (names should be compatible with those in `input`)"
        print "         (the line begin with '#' will be ignored)"
        sys.exit(1)

    opts,args = getopt(argv[1:],'',['sep=','mark='])
    if len(args) != 2:
        print "Error: invalid number of args  ",args
        sys.exit(1)

    sep = ''
    mark_file = ''
    use_names = None
    for opt,val in opts:
        if opt == '--sep':
            sep = val
        elif opt == '--mark':
            mark_file = val
        else:
            print "Error: invalid option ",opt
            sys.exit(1)

    if mark_file != '':
        inf = open(mark_file,'r')
        use_names = []
        for line in inf:
            if line.startswith('#'):
                continue
            use_names.append(line.strip())
        inf.close()

    inf = open(args[0],"r")
    line = inf.readline()
    if sep == "":
        names = line.split()
    else:
        names = line.strip().split(sep)
    names = names[1:]
    if mark_file == '':
        mark = np.ones(len(names),dtype=np.bool)
    else:
        mark = np.zeros(len(names),dtype=np.bool)
        for i in xrange(len(names)):
            if names[i] in use_names:
                mark[i] = True

    outf = open(args[1],"w")
    for line in inf:
        if sep == "":
            line = line.split()
        else:
            line = line.strip().split(sep)
        outf.write(line[0])
        output = np.asarray(line[1:])[mark]
        for i in xrange(output.shape[0]):
            outf.write(" %d:%s"%(i+1,output[i]))
        outf.write("\n")
    inf.close()
    outf.close()

main()

