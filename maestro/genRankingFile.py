import sys
try:
    from schrodinger import structure
except ImportError:
    print r"\n  please run E:\data\virus\maestro\bin\maestro_variables.bat first!"
    sys.exit(1)

def main(argv=sys.argv):
    if len(argv) < 3:
        print "\n    OBJ: to sort ligands according to `r_i_docking_score`"
        print"          resulting ligands' titles will be saved in `out.txt`."
        print "         If several poses have same title, then the one with"
        print "         best `r_i_docking_score` will be kept!"
        print "\n  Usage: %s in1.maegz[...] out.txt"%argv[0]
        print "  in1.maegz[...]: several docked maegz files"
        print "  out.txt: each line will be: title r_i_docking_score\n"
        sys.exit(1)
    
    ligands = {}  #key=title, value=r_i_docking_score
    for infile in argv[1:-1]:
        inf = structure.StructureReader(infile)
        pro = inf.next()
        for mol in inf:
            if ligands.has_key(mol.title):
                if mol.property["r_i_docking_score"] < ligands[mol.title]:
                    ligands[mol.title] = mol.property["r_i_docking_score"]
            else:
                ligands[mol.title] = mol.property["r_i_docking_score"]
    title_score = [(key,value) for key,value in ligands.iteritems()]
    title_score.sort(key=lambda x:x[1])
    outf = open(argv[-1],'w')
    print >>outf, "title r_i_docking_score"
    for title,score in title_score:
        print >>outf, title,score
    outf.close()

main()


