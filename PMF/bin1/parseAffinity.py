'''
#=============================================================================
#     FileName: parseAffinity.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-03-17 15:04:43
#   LastChange: 2014-03-24 16:49:04
#      History:
#=============================================================================
'''
import sys
import re
from math import log10

def main(argv=sys.argv):
    if len(argv) != 3:
        print "\n  Usage: %s input output"%argv[0]
        print "  OBJ: to convert Ki/Kd to log(Ki/Kd, unit: M)\n"
        sys.exit(1)

    pattern = re.compile("\d+\.?\d*")
    #convert = {'M':1, 'uM':10**-3, 'mM':10**-6, 'nM':10**-9, 'pM':10**-12, 'fM':10**-15}
    convert = {"M":0, "mM":-3, "uM":-6, "nM":-9, "pM":-12, "fM":-15}

    inf = open(argv[1],"r")
    outf = open(argv[2],"w")
    line = inf.readline()
    outf.write(line.strip()+"\t-log(Ki/Kd, M)\n")
    count = 0
    for line in inf:
        count += 1
        affinity = line[line.rfind("\t"):].strip()
        result = re.search(pattern,affinity)
        if result is None:
            print "failed to parse line %d"%(count+1)
            outf.write(line)
        else:
            value = float(affinity[result.start():result.end()])
            unit = affinity[result.end():]
            outf.write(line.strip())
            try:
                value = log10(value) + convert[unit]
            except KeyError:
                print "failed to recognize %s in line %d"%(unit, count+1)
                outf.write("\n")
            else:
                outf.write("\t%g\n"%value)
    inf.close()
    outf.close()

main()

