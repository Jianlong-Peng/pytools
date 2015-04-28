
import sys
import os
try:
    from schrodinger import structure
    from schrodinger import structureutil
    from schrodinger.structutils import build
    from schrodinger.infra import mm
    from schrodinger.structutils.minimize import minimize_structure
except ImportError:
    print "\n Please run `E:/data/virus/maestro/bin/maestro_variable.bat` first!!!\n"
    sys.exit(1)

def writeInpFile(outfile,complex_file):
    outf = open(outfile,"w")
    print >>outf, """


[ SET:COMPLEXES]
    VARCLASS Structures
    FILES %s

[ STAGE:PRIME ]
    STAGECLASS              prime.PrimeStage
    INPUTS                  COMPLEXES
    OUTPUTS                 DOCKED_OUT
    PRIME_TYPE              COVALENT_DOCKING
    LIGAND                  X:1
    INCLUDE_RESIDUE         yes
    NUM_OUTPUT_STRUCT       1


[ USEROUTS ]
    USEROUTS DOCKED_OUT
    STRUCTOUT DOCKED_OUT
"""%complex_file
    outf.close()


def calcDistance(xyz1,xyz2):
    distance = 0.
    for i in xrange(3):
        distance += pow(xyz1[i]-xyz2[i],2)
    return distance

def find_closest_atom(atomS, ligand, candidate):
    min_distance = 4   #the closest distance should be less than 4A, around 3A
    min_atom_num = 0
    for i in candidate:
        distance = calcDistance(atomS.xyz, ligand.atom[i].xyz)
        if distance < min_distance:
            min_distance = distance
            min_atom_num = i
    return min_atom_num

def main(argv=sys.argv):
    if len(argv)!=4 and len(argv)!=5:
        print "\n    OBJ: make complex by connecting `S1150` of CYS147 to the closest ligand atom"
        print "         Complexes will be saved in several .maegz files, out_file_1.maegz"
        print "         out_file_2.maegz,... each contains `n` structures"
        print "\n  Usage: makeComplex.py docked_file.maegz out_file.maegz check_file [n]"
        print "  check_file: if given, the closest atom will be one of those in {check_file}"
        print "     from 2nd line on, each line should be `ligand_title reaction_center`"
        print "  [n]: number of complexes each file will contain."
        print "     if not given, all will be written in a single file\n"
        print "  corresponding .inp file(s) will be generated in the same dir as"
        print "  `out_file.maegz`\n"
        sys.exit(1)
    
    assert not os.path.exists(argv[2])
    assert os.path.exists(argv[1])
    
    if len(argv) == 5:
        n = int(argv[4])
        assert n>0
    else:
        n = 0
    
    ligand_atom = {}  #key=title, value=[reaction_atoms,...]
    inf_check = open(argv[3],"r")
    line = inf_check.readline()
    for line in inf_check:
        line = line.split()
        if not ligand_atom.has_key(line[0]):
            ligand_atom[line[0]] = set()
        ligand_atom[line[0]].add(int(line[1]))
    inf_check.close()
    
    basename,ext = os.path.splitext(argv[2])
    
    inf = structure.StructureReader(argv[1])
    receptor = inf.next()
    atomS = receptor.atom[1150]
    count = 0 #number of ligands read
    num = 0   #{num}th complex written to the current output file
    i = 0     #{i}th output file
    out_file_name = ""
    inp_file_name = ""
    if n == 0:
        inp_file_name = "%s.inp"%basename
        out_file_name = argv[2]
    else:
        i = 1
        inp_file_name = "%s_%d.inp"%(basename,i)
        out_file_name = "%s_%d%s"%(basename,i,ext)
    outf = structure.StructureWriter(out_file_name)
    
    for ligand in inf:
        count += 1
        num += 1
        print "> to process ligand `%s`"%(ligand.title)
        
        if n!=0 and num>n:
            outf.close()
            writeInpFile(inp_file_name,out_file_name)
            print "totally %d complexes being written to"%n,out_file_name
            print "corresponding .inp file is:",inp_file_name
            i += 1
            inp_file_name = "%s_%d.inp"%(basename,i)
            out_file_name = "%s_%d%s"%(basename,i,ext)
            outf = structure.StructureWriter(out_file_name)
            num = 1
        #1. find out the closest atom from `S`
        if not ligand_atom.has_key(ligand.title):
            print "Error: ligand `%s` can't be found in `%s`"%(ligand.title,argv[4])
            continue
        reactive_atom_num = find_closest_atom(atomS, ligand, ligand_atom[ligand.title])
        if reactive_atom_num == 0:
            print "Error: none of atom of ligand `%s` is around `S`!!"%ligand.title
            continue
        for atom in ligand.atom:
            atom.chain = 'X'
            atom.resnum = 1
            atom.inscode = ' '
        
        #2. delete hydrogen bound to `S`
        complex_st = receptor.copy()
        receptor_leaving_atom = 0
        for atom in complex_st.atom[1150].bonded_atoms:
            if atom.element == "H":
                receptor_leaving_atom = atom.index
                break
        complex_st.deleteBond(1150, receptor_leaving_atom)
        complex_st.deleteAtoms([receptor_leaving_atom])
        complex_st.extend(ligand)
        complex_st.title = ligand.title
        complex_st.property['s_user_receptor_title'] = receptor.title
        ligand_binding_atom = receptor.atom_total-1 + reactive_atom_num
        #3. delete hydrogen bound to `ligand_binding_atom`
        lig_leaving_atom = 0
        for atom in complex_st.atom[ligand_binding_atom].bonded_atoms:
            if atom.element == "H":
                lig_leaving_atom = atom.index
        complex_st.deleteBond(ligand_binding_atom, lig_leaving_atom)
        complex_st.deleteAtoms([lig_leaving_atom])
        #3. make complex just by connecting `1150` and `ligand_binding_atom`
        complex_st.addBond(ligand_binding_atom, 1150, 1)
        #modify bond length ???????????????????????????
        #the bond length in PDB(code, 3SJO) is 1.80 A
        #complex_st.getBond(ligand_binding_atom, 1150).length = 1.83
        outf.append(complex_st)
    
    if n==0 or num<=n:
        outf.close()
        writeInpFile(inp_file_name,out_file_name)
        print "totally %d complexes being written to"%num,out_file_name
        print "corresponding .inp file is:",inp_file_name
    
main()