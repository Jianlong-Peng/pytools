import sys
import os
try:
    from schrodinger import structure
    from schrodinger import structureutil
except ImportError:
    sys.exit(1)

def get_ligand(complex_st):
    atomS = None
    for atom in complex_st.atom:
        if atom.resnum==147 and atom.element=="S":
            atomS = atom
            break
    assert atomS is not None
    
    lig_atom_break = None
    for atom in atomS.bonded_atoms:
        if atom.chain == "X":
            lig_atom_break = atom
            break
    assert lig_atom_break is not None
    lig_atom_break.property['s_user_breaking_atom'] = '1'
    
    ligand_index = [i for i in xrange(1,complex_st.atom_total+1) if complex_st.atom[i].chain=="X"]
    ligand = complex_st.extract(ligand_index)
    ligand.title = complex_st.title
    atom_break_num = 0
    for atom in ligand.atom:
        if atom.property.has_key('s_user_breaking_atom'):
            atom_break_num = atom.index
    structureutil.add_hydrogens(ligand, atom_list=[atom_break_num])
    
    return ligand



def main(argv=sys.argv):
    if len(argv) != 3:
        print "\n    OBJ: to extract covalent binding ligands"
        print "\n  Usage: %s infile out.maegz"%argv[0]
        print "  infile:"
        print "    .maegz: docked file"
        print "    .txt  : one .maegz file per line"
        sys.exit(1)
    
    assert os.path.exists(argv[1]) and (not os.path.exists(argv[2]))
    
    outf = structure.StructureWriter(argv[2])
    basename,ext = os.path.splitext(argv[1])
    count = 0
    if ext == ".txt":
        inf = open(argv[1],"r")
        for line in inf:
            inf_structures = structure.StructureReader(line.strip())
            for complex_st in inf_structures:
                ligand = get_ligand(complex_st)
                outf.append(ligand)
                count += 1
        inf.close()
    else:
        inf = structure.StructureReader(argv[1])
        for complex_st in inf:
            ligand = get_ligand(complex_st)
            outf.append(ligand)
            count += 1
    outf.close()
    print "Totally extracted %d ligands"%count


main()


