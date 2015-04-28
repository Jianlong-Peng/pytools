'''
#=============================================================================
#     FileName: genSALabel.py
#         Desc: to generate labels used for Yong's structure alert program
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-08-11 10:12:38
#   LastChange: 2012-08-11 10:14:29
#      History:
#=============================================================================
'''

import sys
from math import log10

def main(argv=sys.argv):
    if len(argv) != 3:
        print "Usage: genSALabel.py num_pos num_neg"
        print "to generate labels used for Yong's structure alert program"
        print "generated labels: 0001p, 0010p, 0463n"
        sys.exit(1)

    pos = int(argv[1])
    neg = int(argv[2])

    i = 1
    while i <= pos:
        label = ""
        for j in xrange(4-int(log10(i))-1):
            label += "0"
        label += str(i)
        label += "p"
        print label
        i += 1
    while i<=(pos+neg):
        label = ""
        for j in xrange(4-int(log10(i))-1):
            label += "0"
        label += str(i)
        label += "n"
        print label
        i += 1

main()



