
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

def find_closest_atom(atomS, ligand, match_list):
    min_distance = 1e8
    min_atom_num = 0
    for atom_nums in match_list:
        distance = calcDistance(atomS.xyz, ligand.atom[atom_nums[0]].xyz)
        if distance < min_distance:
            min_distance = distance
            min_atom_num = atom_nums[0]
    return min_atom_num

def main(argv=sys.argv):
    if len(argv) != 4:
        print "\n    OBJ: make complex by connecting `S1150` of CYS147 to the closest atom"
        print "         which matches either `C=CC=O` or `[CH]=O.`"
        print "         Complexes will be saved in several .maegz files, out_file_1.maegz"
        print "         out_file_2.maegz,... each contains `n` structures"
        print "\n  Usage: makeComplex.py docked_file.maegz out_file.maegz n"
        print "  n: number of complexes each file will contain."
        print "     if 0, all will be written to a single file.\n"
        print "  corresponding .inp file(s) will be generated in the same dir as"
        print "  `out_file.maegz`\n"
        sys.exit(1)
    
    assert not os.path.exists(argv[2])
    assert os.path.exists(argv[1])
    
    n = int(argv[3])
    assert n>=0
    
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
        outf = structure.StructureWriter(argv[2])
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
        #1. process alpha,beta-unsaturated ketone. double -> single
        #1.1 find atoms match either `C=CC=O` or `[CH]=O`
        match_list = structureutil.evaluate_smarts(ligand,"C=CC=O")
        match_list.extend(structureutil.evaluate_smarts(ligand,"[CH]=O"))
        if len(match_list) == 0:
            print "  Error: no atom found to match either `C=CC=O` or `[CH]=O`"
            continue
        #1.2 from those matched atoms, selecting the one which is closest to `S` of CYS147
        reactive_atom_num = find_closest_atom(atomS, ligand, match_list)
        #1.3 double -> single
        another_atom_num = 0
        for atom in ligand.atom[reactive_atom_num].bonded_atoms:
            bond = ligand.getBond(reactive_atom_num,atom.index)
            if bond.order == 2:
                another_atom_num = atom.index
                break
        ligand.addBond(reactive_atom_num, another_atom_num, 1)
        structureutil.add_hydrogens(ligand, atom_list=[reactive_atom_num,another_atom_num])
        
        for atom in ligand.atom:
            atom.chain = 'X'
            atom.resnum = 1
            atom.inscode = ' '
        
        #2. make complex just by connecting `1150` and `reactive_atom_num`
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
        complex_st.addBond(ligand_binding_atom, 1150, 1)
        #modify bond length ???????????????????????????
        #the bond length in PDB(code, 3SJO) is 1.80 A
        #complex_st.getBond(ligand_binding_atom, 1150).length = 1.83
        #delete two hydrogens that is connected to `S` and `reactive_atom_num` correspondingly
        lig_leaving_atom = 0
        for atom in complex_st.atom[ligand_binding_atom].bonded_atoms:
            if atom.element == "H":
                lig_leaving_atom = atom.index
        complex_st.deleteBond(ligand_binding_atom, lig_leaving_atom)
        complex_st.deleteAtoms([lig_leaving_atom])
        outf.append(complex_st)
    
    if n == 0:
        outf.close()
    else:
        if num <= n:
            outf.close()
            writeInpFile(inp_file_name,out_file_name)
            print "totally %d complexes being written to"%num,out_file_name
            print "corresponding .inp file is:",inp_file_name
    
main()