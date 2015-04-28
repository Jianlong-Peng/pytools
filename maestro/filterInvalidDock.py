import sys
import os
from math import pow,sqrt
try:
    from schrodinger import structure
except ImportError:
    print "\n  Please run E:\\data\\virus\\maestro\\bin\\maestro_variable.bat first!\n"
    sys.exit(1)

def calcDistance(xyz1, xyz2):
    distance = 0.
    for i in xrange(3):
        distance += pow(xyz1[i]-xyz2[i],2)
    return sqrt(distance)

def find_closest(lig):
    candidate = []
    for atom in lig.atom:
        distance = calcDistance(atom.xyz,(7.38,-12.05,0.06))
        if distance < 1.5:
            candidate.append(atom.index)
    return candidate

def main(argv=sys.argv):
    if len(argv) != 4:
        print "\n  Usage: %s prepare.log docked.maegz out.maegz"%argv[0]
        print "  prepare.log : from 2nd line on, each should be `title reactive_center`"
        print "  docked.maegz: generated by glide docking"
        print "  out.maegz   : only those whose closest atom within (7.38, -12.05, 0.06) r=1A"
        print "                was in `prepare.log` will be saved."
        sys.exit(1)

    assert not os.path.exists(argv[3])
    inf_log = open(argv[1],"r")
    line = inf_log.readline()
    info = {}
    for line in inf_log:
        line = line.split()
        if not info.has_key(line[0]):
            info[line[0]] = set()
        info[line[0]].add(int(line[1]))
    inf_log.close()
    
    inf = structure.StructureReader(argv[2])
    recp = inf.next()
    outf = structure.StructureWriter(argv[3])
    outf.append(recp)
    
    count = 0
    for lig in inf:
        atom_index_list = find_closest(lig)
        find = False
        for atom_index in atom_index_list:
            if atom_index in info[lig.title]:
                find = True
                break
        if find:
            count += 1
            outf.append(lig)
    outf.close()
    print "totally found %d molecules satisfying the criterion `%s`"%(count,argv[1])

main()