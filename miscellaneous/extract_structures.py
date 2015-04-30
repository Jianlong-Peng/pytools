'''
#=============================================================================
#     FileName: extract_structures.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2015-04-09 20:18:19
#   LastChange: 2015-04-09 20:29:39
#      History:
#=============================================================================
'''
import sys
import pybel

def main(argv=sys.argv):
    if len(argv) != 4:
        print "\n  Usage: %s in.sdf in.list out.sdf"%argv[0]
        print ""
        sys.exit(1)

    candidates = [line.strip() for line in open(argv[2],'r')]
    extracted  = [False for _ in candidates]
    outf = pybel.Outputfile("sdf", argv[3])
    inf = pybel.readfile("sdf", argv[1])
    for mol in inf:
        if False not in extracted:
            break
        for i in xrange(len(candidates)):
            if extracted[i]:
                continue
            found = False
            for key,value in mol.data.iteritems():
                if value.lower() == candidates[i].lower():
                    found = True
                    break
            if found:
                extracted[i] = True
                mol.title = candidates[i]
                outf.write(mol)
                break
    inf.close()
    outf.close()

    n = extracted.count(True)
    print "totally extracted %d molecules from %s"%(n, argv[1])
    if n < len(candidates):
        print "The following can't be found"
        for i in xrange(len(extracted)):
            if not extracted[i]:
                print candidates[i]

main()

