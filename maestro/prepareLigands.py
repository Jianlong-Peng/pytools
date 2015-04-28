import sys
import os
try:
    from schrodinger import structure
    from schrodinger import structureutil
    from schrodinger.structutils import build
    from schrodinger.structutils.minimize import minimize_structure
    from schrodinger.infra import mm
except ImportError:
    print r"  Please run E:\data\virus\maestro\bin\maestro_variable.bat first!!"
    sys.exit(1)

purpose = """
  To furthur prepare ligands after applying maestro's `ligprep`
  1. alpha,beta-unsaturated ketone -> ketone
     (C=C double bond -> C-C single bond)
     (all combination(at most 4) of stereoisomers will be generated)
  2. aldehyde -> alcohol
     ([CH]=O -> [CH2]-OH)
  each possible site will be processed and saved as a new ligand.
  corresponding atom index will be saved in a .log file!!!!
"""

def process_ligand(ligand):
    '''
    to process the given ligand
    
    parameter
    =========
    ligand: 
    
    return
    ======
    [(new_ligand,reaction_center),...]
    '''
    result = []  #[(new_ligand,reaction_center),...]
    
    #alpha,beta-unsaturated ketone
    match_list = structureutil.evaluate_smarts(ligand,"C=CC=O")
    for matched_atoms in match_list:
        new_ligand = ligand.copy()
        #double -> single
        new_ligand.addBond(matched_atoms[0],matched_atoms[1],1)
        structureutil.add_hydrogens(new_ligand, atom_list=[matched_atoms[0],matched_atoms[1]])
        minimize_structure(new_ligand)
        result.append((new_ligand,matched_atoms[0]))
        #chirality
        chirality1 = new_ligand.atom[matched_atoms[0]].chirality
        chirality2 = new_ligand.atom[matched_atoms[1]].chirality
        
        #invert chirality of matched_atoms[0]
        if chirality1 in ["R","S"]:
            new_ligand_invert1 = new_ligand.copy()
            #fix `newly added hydrogen which is bonded to matched_atoms[0]`, `alpha-carbon`
            mm.mmct_atom_invert_chirality(new_ligand_invert1, matched_atoms[0], new_ligand.atom_total-1, matched_atoms[1])
            minimize_structure(new_ligand_invert1)
            result.append((new_ligand_invert1,matched_atoms[0]))
        """ Only consider chirality of the reactive atom!!!!!!
        #invert chirality of matched_atoms[1]
        if chirality2 in ["R","S"]:
            new_ligand_invert2 = new_ligand.copy()
            mm.mmct_atom_invert_chirality(new_ligand_invert2, matched_atoms[1], new_ligand.atom_total, matched_atoms[0])
            minimize_structure(new_ligand_invert2)
            result.append((new_ligand_invert2,matched_atoms[0]))
        #invert chirality of both matched_atoms[0] and matched_atoms[1], if possible
        if (chirality1 in ["R","S"]) and (chirality2 in ["R","S"]):
            new_ligand_invert12 = new_ligand.copy()
            mm.mmct_atom_invert_chirality(new_ligand_invert12, matched_atoms[0], new_ligand.atom_total-1, matched_atoms[1])
            mm.mmct_atom_invert_chirality(new_ligand_invert12, matched_atoms[1], new_ligand.atom_total, matched_atoms[0])
            minimize_structure(new_ligand_invert12)
            result.append((new_ligand_invert12,matched_atoms[0]))        
        """
    #aldehyde
    match_list = structureutil.evaluate_smarts(ligand,"[CH]=O")
    for matched_atoms in match_list:
        new_ligand = ligand.copy()
        #double -> single
        new_ligand.addBond(matched_atoms[0],matched_atoms[1],1)
        structureutil.add_hydrogens(new_ligand, atom_list=[matched_atoms[0],matched_atoms[1]])
        minimize_structure(new_ligand)
        result.append((new_ligand,matched_atoms[0]))
        #no chirality!
    
    return result


def main(argv=sys.argv):
    global purpose
    if len(argv) != 4:
        print "\n  Description:"
        print purpose
        print "\n  Usage: %s in.maegz out.maegz out.log\n"%argv[0]
        print "  in.maegz : ligands prepared via `ligprep`"
        print "  out.maegz: to save newly generated ligands"
        print "  out.log  : to save atom number of each ligand's reaction center\n"
        sys.exit(1)
    
    assert os.path.exists(argv[1])
    assert not os.path.exists(argv[2])
    
    inf = structure.StructureReader(argv[1])
    outf = structure.StructureWriter(argv[2])
    outf_log = open(argv[3],"w")
    print >>outf_log, "ligand_title reaction_center"
    count_in = 0
    count_out = 0
    for ligand in inf:
        count_in += 1
        new_ligands = process_ligand(ligand)
        if len(new_ligands) == 0:
            print "ligand `%s` doesn't match either `C=CC=O` or `[CH]=O`"%(ligand.title)
            continue
        for st,atom_index in new_ligands:
            count_out += 1
            outf.append(st)
            print >>outf_log, st.title,atom_index
    outf.close()
    outf_log.close()
    print "Totally read %d ligands. After processing, generated %d ligands!"%(count_in,count_out)

main()
