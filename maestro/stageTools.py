'''
#=============================================================================
#     FileName: stageTools.py
#         Desc: make complex from protein(Cys) and ligand(alpha,beta-unsaturated ketone)
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-06-21 14:15:25
#   LastChange: 2013-06-21 15:51:35
#      History:
#=============================================================================
'''
import sys
from schrodinger import structureutil
from schrodinger.structutils import build
from schrodinger.structutils.minimize import minimize_structure
from schrodinger.infra import mm
import schrodinger.infra.mmbitset as mmbitset


def deleteReceptorLeavingGroup(receptor, leaving_atom_num, staying_atom_num):
    '''
    '''
    st = receptor.copy()
    if st.atom[leaving_atom_num].element == 'H':
        return st
    print "  Error<deleteReceptorLeavingGroup>: currently leaving atom should be `H`, but %s given!"%(st.atom[leaving_atom_num].element)
    return None


def modifyLigandDoubleBond(ligand, reactive_atom_num):
    '''
    to convert double bond to single one!

    parameters
    ==========
    ligand: schrodinger.structure.Structure
    reactive_atom_num: int. atom which will be bound to the protein

    returns
    =======
    list of modified ligand
    '''
    st = ligand.copy()
    reactive_atom = st.atom[reactive_atom_num]

    orig_nbors = [atom for atom in reactive_atom.bonded_atoms]
    atom_pair_in_ring = []
    for i in xrange(len(orig_nbors)-1):
        for j in xrange(i+1, len(orig_nbors)):
            if st.inRing(orig_nbors[i], orig_nbors[j]):
                atom_pair_in_ring.append((int(orig_nbors[i]),int(orig_nbors[j])))
    
    #get the atom which is connected to `reactive_atom` by double bond!
    another_atom_num = 0
    #nonH_nbor_atom_num = 0  #any nonhydrogen or nonhydrogen with less fragments attached ????
    for atom in orig_nbors:
        bond = st.getBond(reactive_atom_num,atom.index)
        if bond.order == 2:
            another_atom_num = atom.index
            break
    
    if another_atom_num == 0:
        print "  Error<modifyLigandDoubleBond>: no double bond found for the given atom",reactive_atom_num
        return None
            
    #modify bond type from `double` to `single`
    st.addBond(reactive_atom_num,another_atom_num,1)
    #add Hydrogen both to leaving and staying atoms
    structureutil.add_hydrogens(st, atom_list=[reactive_atom_num,another_atom_num])
    #optimize the coordinates   2013-06-24
    minimize_structure(st)
    
    ligands = [st.copy()]
    
    #modify  2013-06-24  21:45
    if atom_pair_in_ring == []:
        new_h_atom = [atom for atom in reactive_atom.bonded_atoms if atom not in orig_nbors]
        assert len(new_h_atom) == 1
        fixed_atom1 = new_h_atom[0]
        for nbor in orig_nbors:
            fixed_atom2 = int(nbor)
            tmp_st = st.copy()
            mm.mmct_atom_invert_chirality(tmp_st, reactive_atom_num, fixed_atom1, fixed_atom2)
            minimize_structure(tmp_st)
            ligands.append(tmp_st)
    else:
        #if there exists ring, fix those neighbor atoms who are in the ring!!!!
        fixed_atom1,fixed_atom2 = atom_pair_in_ring[0]
        tmp_st = st.copy()
        mm.mmct_atom_invert_chirality(tmp_stm, reactive_atom_num, fixed_atom1, fixed_atom2)
        minimize_structure(tmp_st)
        ligands.append(tmp_st)

    return ligands


def makeComplex(receptor, pro_leaving_atom, pro_staying_atom, ligand, lig_reactive_atom):
    '''
    parameters
    ==========
    receptor, ligand: schrodinger.structure.Structure
    lig_reactive_atom: int. which will be bound to the `pro_staying_atom`
    '''
    receptor_st = deleteReceptorLeavingGroup(receptor, pro_leaving_atom, pro_staying_atom)
    if receptor_st is None:
        return None
    
    ligands_st = modifyLigandDoubleBond(ligand, lig_reactive_atom)
    if ligands_st is None:
        return None
    
    lig_reactive_atom += receptor_st.atom_total
    complexes = []
    for ligand_st in ligands_st:
        #mark ligand
        #so LIGAND of *.inp must be "X:1"
        for atom in ligand_st.atom:
            atom.chain = 'X'
            atom.resnum = 1
            atom.inscode = ' '

        complex_st = receptor_st.copy()
        complex_st.extend(ligand_st)

        complex_st.title = ligand_st.title
        complex_st.property['s_user_receptor_title'] = receptor_st.title
    
        try:
            complex_st.write('complex_st.mae')
            #if there exist several possible connection ??????????
            renumber_dict = build.connect(complex_st, [pro_staying_atom], [lig_reactive_atom])
        except mm.MmException:
            print "  Error<makeComplex>: failed to connect receptor atom %d to ligand atom %d"%(pro_staying_atom, lig_reactive_atom)
            sys.exit(1)
        else:
            complexes.append(complex_st)
    
    return complexes


"""
def makeComplex(receptor, pro_leaving_atom, pro_staying_atom, ligand, lig_reactive_atom):
    '''
    parameters
    ==========
    receptor, ligand: schrodinger.structure.Structure
    lig_reactive_atom: int. which will be bound to the `pro_staying_atom`
    '''
    receptor_st = deleteReceptorLeavingGroup(receptor, pro_leaving_atom, pro_staying_atom)
    if receptor_st is None:
        return None

    ligand_st = modifyLigandDoubleBond(ligand, lig_reactive_atom)
    if ligand_st is None:
        return None
    #mark ligand
    #so LIGAND of *.inp must be "X:1"
    for atom in ligand_st.atom:
        atom.chain = 'X'
        atom.resnum = 1
        atom.inscode = ' '

    complex_st = receptor_st.copy()
    complex_st.extend(ligand_st)
    lig_reactive_atom += receptor_st.atom_total

    complex_st.title = ligand_st.title
    complex_st.property['s_user_receptor_title'] = receptor_st.title
    
    try:
        complex_st.write('complex_st.mae')
        #if there exist several possible connection ??????????
        renumber_dict = build.connect(complex_st, [pro_staying_atom], [lig_reactive_atom])
    except mm.MmException:
        print "  Error<makeComplex>: failed to connect receptor atom %d to ligand atom %d"%(pro_staying_atom, lig_reactive_atom)
        return None
    else:
        #modify. 2013-06-24 20:15
        lig_atom = complex_st.atom[renumber_dict[lig_reactive_atom]]
        rec_atom = complex_st.atom[renumber_dict[pro_staying_atom]]
        complexes = [complex_st.copy()]
        chirality = lig_atom.chirality
        #invert chirality of `lig_reactive_atom`
        if chirality in ['R','S']:
            lig_atom.label_user_text = chirality
            lig_atom.label_format = "%UT"
            lig_atom.label_color = 13
            
            fixed_atom1 = int(rec_atom)
            fixed_atom2 = None
            
            ligand_atoms = [atom for atom in lig_atom.bonded_atoms if atom != rec_atom]
            
            for atom in ligand_atoms:
                in_ring = complex_st.inRing(rec_atom,atom)
                if in_ring:
                    fixed_atom2 = int(atom)
                    break
            if not fixed_atom2:
                fixed_atom2 = int(ligand_atoms[0])
            
            mm.mmct_atom_invert_chirality(complex_st, int(lig_atom), fixed_atom1, fixed_atom2)
            lig_atom.label_user_text = lig_atom.chirality
            complexes.append(complex_st.copy())
        return complexes
"""

"""
def deleteLigandLeavingGroup(st, leaving_atom_num):
        '''
        Called for both ligand structures.
        Finds the leaving group atoms and deletes them.
        Returns modified structure and atom number which
        the leaving group was bound to.
        '''
        
        st = st.copy() # Do not modify the original

        leaving_atom = st.atom[leaving_atom_num]
        
        leaving_group_atoms = []
        
        if leaving_atom.element == 'H':
            staying_atom = leaving_atom.bond[1].atom2
        elif len(leaving_atom.bond) == 1:
            staying_atom = leaving_atom.bond[1].atom2
            # Leaving atom will get replaced with a hydrogen:
            leaving_group_atoms.append( int(leaving_atom) )
        
        else:
            # Leaving group is not a single atom. Find other atoms:
            # Break the single/double bond, delete the leaving group,
            # and add hydrogen in place of the leaving group.
            
            # Follow out each branch from the leaving atom.
            # Staying atom is the one whose branch does not terminate:
            
            most_atom_n = None
            most_atom_n_atoms = 0
            
            for n in leaving_atom.bonded_atoms:
                # Determine number of atoms in st that will move along with
                # atom n when leaving_atom is fixed. 
                bs = mmbitset.Bitset(size=st.atom_total)
                mm.mmct_atom_get_moving(st, leaving_atom, st, n, bs)
                num_branch_atoms = bs.count()
                
                if num_branch_atoms > most_atom_n_atoms:
                    most_atom_n = int(n)
                    most_atom_n_atoms = num_branch_atoms
            
            staying_atom = st.atom[most_atom_n]

            st.deleteBond(int(leaving_atom), int(staying_atom))
            leaving_group_atoms = map(int, st.getMoleculeAtoms(leaving_atom))
        
        
        # Delete leaving group atoms and renumber the staying atom:
        if leaving_group_atoms:
            if int(staying_atom) in leaving_group_atoms:
                st.write('tmpligand.mae')
                msg = "ERROR: Atom %i in ligand can not be a leaving group. Please change the SMARTS pattern and SMARTS position so that matching atom is the atom in the leaving group that is directly bound to the rest of the ligand." % leaving_atom.index
                msg += "\n  Leaving atom: %s" % get_atom_str(leaving_atom)
                msg += "\n  Staying atom: %s" % get_atom_str(staying_atom)
                msg += "Structure was written to: tmpligand.mae"
            
            staying_atom.property['i_tmp_staying_atom'] = 1
            st.deleteAtoms(leaving_group_atoms)
            result = structureutil.evaluate_asl(st, "atom.i_tmp_staying_atom 1")
            staying_atom_num = result[0]
            del st.atom[staying_atom_num].property['i_tmp_staying_atom']
            
            # Add a hydrogen to the staying atom:
            structureutil.add_hydrogens(st, atom_list=[staying_atom_num])
        else:
            staying_atom_num = int(staying_atom)
        
        
        return st, staying_atom_num


def deleteReceptorLeavingGroup(st, leaving_atom_num, staying_atom_num):
        '''
        Called for both receptor structure.
        Finds the leaving group atoms and deletes them.
        Returns modified structure and the renumbered staying atom.
        '''
        
        st = st.copy() # Do not modify the original

        staying_atom = st.atom[staying_atom_num]
        leaving_atom = st.atom[leaving_atom_num]
        
        
        leaving_group_atoms = []

        if leaving_atom.element == 'H':
            pass # the hydrogen will later get replaced
        elif len(leaving_atom.bond) == 1:
            # Leaving atom will get replaced with a hydrogen:
            leaving_group_atoms.append( int(leaving_atom) )
        
        else:
            # Leaving group is not a single atom. Find other atoms:
            # Break the single/double bond. The molecule created on the
            # leaving side is the leaving group. Remove it and add a
            # hydrogen in place of it.
            
            st.deleteBond(int(leaving_atom), int(staying_atom))
            leaving_group_atoms = map(int, st.getMoleculeAtoms(leaving_atom))
                
        
        # Delete leaving group atoms and renumber the staying atom:
        if leaving_group_atoms:
            if staying_atom.index in leaving_group_atoms:
                st.write('tmpreceptor.mae')
                return None,0
            
            staying_atom.property['i_tmp_staying_atom'] = 1
            st.deleteAtoms(leaving_group_atoms)
            result = structureutil.evaluate_asl(st, "atom.i_tmp_staying_atom 1")
            staying_atom_num = result[0]
            del st.atom[staying_atom_num].property['i_tmp_staying_atom']
            
            # Add a hydrogen to the staying atom:
            structureutil.add_hydrogens(st, atom_list=[staying_atom_num])
        else:
            staying_atom_num = int(staying_atom)
        
        
        return st, staying_atom_num


def makeComplex(receptor,receptor_leaving_atom,receptor_staying_atom,ligand,ligand_reactive_atom):
        
        # Receptor is already pre-processed (leaving group deleted).
        receptor_st,staying_rec_atom_num = deleteReceptorLeavingGroup(receptor, receptor_leaving_atom, receptor_staying_atom)
        if receptor_st is None:
            return None
        
        # Prepare the ligand (delete leaving group):
        # Atom of the ligand which will leave:
        leaving_lig_atom_num = 0
        for atom in ligand.atom[ligand_reactive_atom].bonded_atoms:
            bond = ligand.getBond(ligand_reactive_atom,atom.index)
            if bond.order == 2:
                leaving_lig_atom_num = atom.index
                break
        assert leaving_lig_atom_num > 0
        ligand_st, staying_lig_atom_num = deleteLigandLeavingGroup(ligand, leaving_lig_atom_num)
        #mark ligand
        #so LIGAND of *.inp must be "X:1"
        for atom in ligand_st.atom:
            atom.chain = 'X'
            atom.resnum = 1
            atom.inscode = ' '

        # Make the complex:
        complex_st = receptor_st.copy()
        complex_st.extend(ligand_st)
        staying_lig_atom_num += receptor_st.atom_total
        
        # Ev:74674 Take title from ligand, and put receptor title into prop:
        complex_st.title = ligand_st.title
        complex_st.property["s_user_receptor_title"] = receptor_st.title
        
        
        try:
            complex_st.write('complex_st.mae')
            renumber_dict = build.connect(complex_st, [staying_rec_atom_num], [staying_lig_atom_num])
        except mm.MmException:
            complex_st.write('tmpcomplex.mae')
            exit(1)
                
        # Ligand and receptor atoms in the complex:
        lig_atom = complex_st.atom[ renumber_dict[staying_lig_atom_num] ]
        rec_atom = complex_st.atom[ renumber_dict[staying_rec_atom_num] ]
        
        # FIXME: Implement Ev:62299 hybridizations
        
        return complex_st
"""
