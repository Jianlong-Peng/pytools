
import sys
try:
    from schrodinger import structure
except ImportError:
    print "\n  Please run E:\\data\\virus\\maestro\\bin\\maestro_variable.bat first!!"
    sys.exit(1)

def main(argv=sys.argv):
    if len(argv)!=2 and len(argv)!=3:
        print "\n  Usage: %s in.maegz [mol_title]"%argv[0]
        print "  in.maegz: file generated via glide dock"
        print "  [mol_title]: if given, only those molecules are extracted!\n"
        print "  title of each ligand will be used as output file name\n"
        sys.exit(1)
    
    if len(argv) == 4:
        mol_title = argv[3]
    else:
        mol_title = None
    
    inf = structure.StructureReader(argv[1])
    receptor = inf.next()
    num = 0
    for ligand in inf:
        write_ligand = False
        if (mol_title is None) or (ligand.title == mol_title):
            num += 1
            write_ligand = True
        if write_ligand:
            complex_st = receptor.copy()
            complex_st.extend(ligand)
            title = ligand.title.replace("/","_")
            complex_st.write("%s.pdb"%title)
    print "Totally %d pdb files were generated!"%num


main()

