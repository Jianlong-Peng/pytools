'''
#=============================================================================
#     FileName: tools.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-06-03 16:13:50
#   LastChange: 2014-06-26 14:57:11
#      History:
#=============================================================================
'''
from math import sqrt
import os

amino_acids_standard = set(['GLY','ALA','VAL','LEU','ILE','PHE','TRP','TYR','ASP',\
        'HIS','ASN','GLU','LYS','GLN','MET','ARG','SER','THR','CYS','PRO'])
nucleotide_standard = set(["DC","DG","DA","DU","DT","DI","C","G","A","U","I"])
#the following atomic mass is from 
#  http://en.wikipedia.org/wiki/Relative_atomic_mass
atomic_mass = {"H":1.008, "C":12.01, "N":14.01, "O":16.00, "F":19.00,\
        "NA":22.99,"MG":24.31,"AL":26.98,"SI":28.09,"P":30.97,"S":32.06,\
        "CL":35.45,"K":39.10,"CA":40.08,"MN":54.94,"FE":55.85,"CU":63.55,\
        "ZN":65.38,"BR":79.90,"I":126.90,"AG":107.87,"PT":195.08,"AU":196.97,\
        "HG":200.59}

common_elements = set(['C','N','O','F','P','S','Cl','Br','I','H'])
min_num_hvy_atoms = 6
max_mw = 1000.
max_peptide_len = 10
max_nucleotide_len = 3
clash_dist = 2.0

#the following 'sepcial' heterogen molecules were extracted from
#  Wang R.X. et al. J. Med. Chem. 2005, 48:4111-4119
#complex contains any of the heterogens will be rejected
special_molecules = set([\
        #coenzyme A family
        'COA','2CP','3CP','3HC','4CA','4CO','ACO','AMX',\
        'BCA','CA3','CA5','CAA','CAO','CIC','CMC','CMX',\
        'CO8','COD','COF','COS','COT','CS8','DAK','DCA',\
        'DCC','FAM','FCX','HAX','HMG','HXC','LYX','MCA',\
        'MCD','MDE','MLC','MYA','NHM','NMX','SCA','SCD',\
        'SCO'," UQ",\
        #Heme (metal-containing protoporphyrin) Family
        'HEM','1CP','B12','BCB','BCL','BLA','BLV','BPB',\
        'BPH','CCH','CL1','CL2','CLA','CLN','CNC','COB',\
        'COH','CON','COJ','COY','CP3','DDH','DEU','DHE',\
        'F43','FEC','HAS','HDD','HDM','HE6','HEA','HEB',\
        'HEC','HEG','HEO','HES','HEV','HIF','HNI','MHM',\
        'MMP','MP1','PC3','PCU','PNI','POR','PP9','SRM',\
        'ZEM','ZHN',\
        #NAD (nicotinamide adenine Dinucleotide) Family
        'NAD','ADJ','CAN','CND','DND','NAC','NAE','NAH',\
        'NAI','NAJ','NAP','NAQ','NAX','NBD','NBP','NDA',\
        'NDC','NDO','NDP','NHD','NHO','ODP','PAD','SND',\
        'TAP','ZID',\
        #FAD (flavin adenine dinucleotide) family
        'FAD','6FA','FAA','FAB','FAE','FAS','FDA','FMA',\
        'FMN','FNS','MGD','RFL'])
#keep??!!
nucleotides = set([\
        'AMP','ADP','ANP','ATP','CMP','CDP','CTP','GMP',\
        'GDP','GNP','GTP','2GP','TMP','TDP','TTP','UMP',\
        'UDP','UTP','PSU'])
junk_molecules = set([\
        #"Junk" molecules
        'BMA','BOG','C8E','CIT','CRY','DTT','EPE','F6P',\
        'FUC','GAL','GLC','GOL','HED','LDA','LI1','MAL',\
        'MAN','MES','MPD','MYR','NAG','NGA','PEG','PG4',\
        'POP','PYR','SPM','TRS','XYS'])

# solvent_inorganic_molecules and ion_molecules are invalid ligands
# from analysis of PDB database
solvent_inorganic_molecules = set([\
        "SO4","EDO","DMS","FMT","ACY","BME","NO3","SX",\
        "SF4","ACE","PGE","IPA","FES","EOH","SCN","CMP",\
        "FLC","CO3","NH2","MLI","OXY","DMF","P6G","TLA",\
        "CAC","URE","AZI","HEZ","PGO","NH4","MOH","CYN",\
        "F3S","BCT","MLA","NO","OH","NO2","O","PI","FOR"\
        "XCC","PO4"])

ion_molecules = set([\
        "CL","ZN","CA","MG","NA","MN","IOD","K","CD","FE",\
        "CU","BR","NI","CO","FE2","HG","XE","CU1","CS",\
        "PT","BA","AU","GD","YB","SR","PB","RB","TL","U1",\
        "LI","F","SM","W","PD","LU","Y1","MO","TB","HO",\
        "EU","KR","DY","RU","AG","SE","LA","OS","AL","AR",\
        "TE","SB","RH","CR","CE","GA","IR"])


class Atom:
    def __init__(self, atom_no, atom_name, element, x, y, z):
        self.atom_no = atom_no
        self.atom_name = atom_name
        #slef.res_name = None
        #self.chain = None
        #self.res_no = None
        self.coords = [x,y,z]
        self.element = element
        self.radius = None

class Residue:
    def __init__(self, res_name, chain, res_no):
        self.atoms = []
        self.res_name = res_name
        self.chain = chain
        self.res_no = res_no

    def get_name(self):
        return self.res_name

    def add_atom(self, atom):
        self.atoms.append(atom)

    def chain_id(self):
        return self.chain

    def __iter__(self):
        self._count = -1
        return self

    def next(self):
        self._count += 1
        if self._count < len(self.atoms):
            return self.atoms[self._count]
        else:
            raise StopIteration


#for polypeptides or polynucleotide
class Polymer:
    def __init__(self):
        self.residues = []   #list of <type Residue>

    def get_name(self):
        return "_".join([item.res_name for item in self.residues])

    def add_residue(self, res):
        self.residues.append(res)

    def is_polypeptide(self):
        for res in self.residues:
            if res.res_name in amino_acids_standard:
                return True
        return False

    def is_polynucleotide(self):
        for res in self.residues:
            if res.res_name in nucleotide_standard:
                return True
        return False

    def chain_id(self):
        if len(self.residues) == 0:
            return ""
        else:
            return self.residues[0].chain

    def __iter__(self):
        self.i = 0
        self.j = -1
        return self

    def next(self):
        self.j += 1
        while self.i<len(self.residues) and self.j>=len(self.residues[self.i].atoms):
            self.i += 1
            self.j = 0
        if self.i<len(self.residues) and self.j<len(self.residues[self.i].atoms):
            return self.residues[self.i].atoms[self.j]
        else:
            raise StopIteration



# to extract contents from given pdb file
#     HEADER, COMPND, TITLE, ATOM, TER, HETATM, CONECT
# input : inf - <type 'file'>
# return: pro_chain, ligands, waters, conect - list of <type string>, main contents
#         conect_dict - dict, key=string, value=[string,...]
def read_pdb_file(inf):
    pro_chain = []
    ligands = []
    waters = []
    conect = []
    conect_dict = {}
    line_no = 0
    inf.seek(0)
    #1. read ATOM/HETATM blocks
    line = inf.readline()
    line_no += 1
    while line!="" and not (line.startswith("ATOM") or line.startswith("HETATM")):
        line = inf.readline()
        line_no += 1
    if line == "":
        print "    Error: can't find 'ATOM' or 'HETATM'"
        return [],[],[],[],{}
    while line != "":
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            break
        chain = []
        res = [[]]
        index_list = [" "]
        occupancy = [0]
        #res = [[] for _ in xrange(5)]
        #index_list = [" ", "A", "B", "C", "D"]
        #occupancy = [0 for _ in xrange(5)]
        prev_res_num = -10000
        while line != "":
            if line.startswith("ANISOU"):
                line = inf.readline()
                line_no += 1
                continue
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                break
            res_num = int(line[22:26])
            if res_num != prev_res_num:
                #print res_num
                if sum(occupancy) != 0:
                    if len(occupancy)==1 or max(occupancy[1:])==0.:
                        i = 0
                    else:
                        i = occupancy[1:].index(max(occupancy[1:])) + 1
                    chain.extend(res[i])
                    res = [[]]
                    index_list = [" "]
                    occupancy = [0]
                    #res = [[] for _ in xrange(5)]
                    #occupancy = [0 for _ in xrange(5)]
                prev_res_num = res_num
            if line[16] not in index_list:
                res.append([item for item in res[0]])
                index_list.append(line[16])
                occupancy.append(0)
            index = index_list.index(line[16])
            line = line[:16]+" "+line[17:]    #delete 'alternate location indicator' ???????
            if index == 0:
                map(lambda item: item.append(line), res)
            else:
                res[index].append(line)
            occupancy[index] = float(line[54:60])
            line = inf.readline()
            line_no += 1
        #suppose to be protein chain without TER
        #if line == "":
        #    print "  Error: there is something wrong in ATOM/HETATM block"
        #    return [],[],[],[],{}
        if sum(occupancy) == 0.:
            print "    Error: there is something wrong with 'occupancy'. L%d"%line_no
            return [],[],[],[],{}
        if len(occupancy)==1 or max(occupancy[1:])==0.:
            i = 0
        else:
            i = occupancy[1:].index(max(occupancy[1:])) + 1
        chain.extend(res[i])
        #protein
        if line=="" or line.startswith("TER"):
            for item in chain:
                if item[17:20].strip() in nucleotide_standard:
                    print "    Error: contains DNA or RNA"
                    return [],[],[],[],{}
                pro_chain.append(item)
            if line != "":
                pro_chain.append(line)  #TER...
        #heterogens
        else:
            for item in chain:
                if item[17:20] == "HOH":
                    waters.append(item)
                else:
                    ligands.append(item)
        #next block
        if line.startswith("TER"):
            line = inf.readline()
            line_no += 1
    #4. try to read CONECT records
    while line!="" and not line.startswith("CONECT"):
        line = inf.readline()
        line_no += 1
    #no CONECT records
    if line == "":
        return pro_chain,ligands,waters,conect,conect_dict
    #read CONECT records - fill 'conect' and 'conect_dict'
    while line!="" and line.startswith("CONECT"):
        i = int(line[6:11])
        i1 = line[11:16].strip()
        if i1 == "":
            line = inf.readline()
            line_no += 1
            continue
        i2 = line[16:21].strip()
        i3 = line[21:26].strip()
        i4 = line[26:31].strip()
        conect_dict[i] = [int(i1)]
        if i2 != "":
            conect_dict[i].append(int(i2))
        if i3 != "":
            conect_dict[i].append(int(i3))
        if i4 != "":
            conect_dict[i].append(int(i4))
        conect.append(line)
        line = inf.readline()
        line_no += 1

    return pro_chain,ligands,waters,conect,conect_dict


#called by 'process'
#OBJ: if contains any nucleotide (standard)
#need a more quick way !!!!!!!!!!!!!!!!!!
def check_contain_nucleotide(pro_chain):
    for line in pro:
        if line[17:20].strip() in nucleotide_standard:
            print "    complex contains DNA/RNA ",line[17:20]
            return False
    return True

#called by 'construct_ligands'
#OBJ: to check if the two ligand is covalently bound
#input:  ligand1, ligand2    <type Residue>
#        conect_dict         dict, key=id, value=[id,...]
#return: True or False
def check_ligand_covalent_bound(ligand1,ligand2, conect_dict):
    atom_no_1 = [atom.atom_no for atom in ligand1]
    atom_no_2 = set([atom.atom_no for atom in ligand2])
    for no in atom_no_1:
        if not conect_dict.has_key(no):
            continue
        for conect_no in conect_dict[no]:
            if conect_no in atom_no_2:
                return True
    return False


#called by 'process'
#input:   ligands  list of <type string>
#         conect_dict    dict, key=id, value=[id,...]
#return:  final_ligands    [<type Residue>, <type Polymer>, ...]
#use chain_identifier(L22) and residue_number(L23-26) to determine
#a ligand, followed by combination if any of them were covalently bound.
def construct_ligands(ligands, conect_dict):
    #1. extract each residue(ligand)
    all_ligands = []
    each_lig = None
    prev_chain = ""
    prev_res_no = -10000
    for line in ligands:
        atom_no = int(line[6:11])
        atom_name = line[12:16].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        element = line[76:78].strip()
        each_atom = Atom(atom_no, atom_name, element,x,y,z)
        res_name = line[17:20].strip()
        chain = line[21]
        res_no = int(line[22:26])
        if chain!=prev_chain or res_no!=prev_res_no:
            if each_lig is not None:
                all_ligands.append(each_lig)
            each_lig = Residue(res_name, chain, res_no)
            prev_chain = chain
            prev_res_no = res_no
        each_lig.add_atom(each_atom)
    if each_lig is not None:
        all_ligands.append(each_lig)
    #2. re-check all the residues, and merge those covalently bound
    final_ligands = []
    combined = [False for _ in all_ligands]
    for i in xrange(len(all_ligands)):
        if combined[i]:
            continue
        combine_list = [i]
        for j in xrange(i+1, len(all_ligands)):
            if combined[j]:
                continue
            if check_ligand_covalent_bound(all_ligands[i],all_ligands[j], conect_dict):
                combine_list.append(j)
                combined[j] = True
        if len(combine_list) == 1:
            final_ligands.append(all_ligands[i])
        else:
            each = Polymer()
            for k in combine_list:
                each.add_residue(all_ligands[k])
            final_ligands.append(each)

    return final_ligands


#called by 'check_valid_ligand'
#input: heterogens     [<type Residue>, <type Polymer>, ...]
#       ligand_valid   [0,1,2,...]
#return:
#  if contains special_molecules, return True
#  otherwise, return false
#junk_molecule will be checked as invalid
def check_special(heterogens, ligand_valid):
    for i in xrange(len(heterogens)):
        if ligand_valid[i] != 1:
            continue
        heterogen = heterogens[i]
        if isinstance(heterogen, Residue):
            if heterogen.res_name in special_molecules:
                print "    PDB contains special ligand: ",heterogen.res_name
                return True
            if heterogen.res_name in junk_molecules:
                print "    PDB contains 'junk' ligand: ",heterogen.res_name
                ligand_valid[i] = 0
        else:
            num_junk = 0
            for res in heterogen.residues:
                if res.res_name in special_molecules:
                    print "    PDB contains special ligand: ",res.res_name
                    return True
                if res.res_name in junk_molecules:
                    num_junk += 1
            if num_junk == len(heterogen.residues):
                print "    PDB contains 'junk' ligand: ",[res.res_name for res in heterogen.residues]
                ligand_valid[i] = 0
    return False


#called by 'check_valid_ligand'
#input: heterogens   [<type Residue>, <type Polymer>, ...]
#       ligand_valid [0,1,2,...]   (0 - invalid; 1 - valid; 2 - ions)
def check_ligand_solvent_ion(heterogens, ligand_valid):
    for i in xrange(len(heterogens)):
        if isinstance(heterogens[i], Residue):
            if heterogens[i].res_name in solvent_inorganic_molecules:
                print "    contains solvent/inorganic: ",heterogens[i].res_name
                ligand_valid[i] = 0
            elif heterogens[i].res_name in ion_molecules:
                print "    contains metal/ion: ",heterogens[i].res_name
                ligand_valid[i] = 2
            else:
                pass


#called by 'check_valid_ligand'
#OBJ  :  to check if contains uncommon elements or MW<=1000
#        or < 6 non-hydrogen atoms
#input:  heterogens     [<type Residue>, <type Polymer>, ...]
#        ligand_valid   [0,1,2,...]   (0 - invalid; 1 - valid; 2 - ions)
def check_common_element(heterogens, ligand_valid):
    for i in xrange(len(heterogens)):
        if ligand_valid[i] != 1:
            continue
        MW = 0
        num_hvy = 0
        for atom in heterogens[i]:
            if atom.element not in common_elements:
                ligand_valid[i] = 0
                print "    ligand(%s) contain uncommon element: %s"%(heterogens[i].get_name(), atom.element)
                break
            MW += atomic_mass[atom.element]
            if atom.element != "H":
                num_hvy += 1
        if ligand_valid[i] and (num_hvy<min_num_hvy_atoms or MW>max_mw):
            ligand_valid[i] = 0
            print "    ligand(%s): number of heavy atoms %d, MW %g"%(heterogens[i].get_name(), num_hvy, MW)


#called by 'check_valid_ligand'
#OBJ: 
#      - for polypeptide, <= 10 residues
#      - for polynucleotide, <= 3 residues
#input:  heterogens     [<type Residue>, <type Polymer>, ...]
#        ligand_valid   [0,1,2,...]   (0 - invalid; 1 - valid; 2 - ions)
def check_polymer(heterogens, ligand_valid):
    for i in xrange(len(heterogens)):
        if ligand_valid[i]!=1 or (not isinstance(heterogens[i], Polymer)):
            continue
        if heterogens[i].is_polypeptide() and len(heterogens[i].residues)>max_peptide_len:
            ligand_valid[i] = 0
            print "    complex contains polypetide with %d residues"%len(heterogens[i].residues)
        elif heterogens[i].is_polynucleotide() and len(heterogens[i].residues)>max_nucleotide_len:
            ligand_valid[i] = 0
            print "    complex contains polynucleotide with %d residues"%len(heterogens[i].residues)
        else:
            pass


#called by 'process'
#input:   heterogens   [<type Residue>, <type Polymer>, ...]
def check_valid_ligand(heterogens):
    ligand_valid = [1 for i in xrange(len(heterogens))]  #[0,1,2,...]  (0 - invalid; 1 - valid; 2 - ions)
    #if contains special heterogens
    if check_special(heterogens, ligand_valid):
        print "    no more further checking"
        return [0 for i in xrange(len(heterogens))]
    if 1 not in ligand_valid:
        return [0 for i in xrange(len(heterogens))]
    #if contains solvents/inorganic or ions
    check_ligand_solvent_ion(heterogens, ligand_valid)
    if 1 not in ligand_valid:
        return [0 for i in xrange(len(heterogens))]
    #if contains only common elements and MW<=1000
    #and more than 6 non-hydrogen atoms
    check_common_element(heterogens, ligand_valid)
    if 1 not in ligand_valid:
        return [0 for i in xrange(len(heterogens))]
    #check polymer 
    # - for polypeptide, <= 10 residues
    # - for polynucleotide, <= 3 residues
    check_polymer(heterogens, ligand_valid)
    if 1 not in ligand_valid:
        return [0 for i in xrange(len(heterogens))]

    return ligand_valid


calc_distance = lambda coord1, coord2: \
        sqrt(pow(coord1[0]-coord2[0], 2) + pow(coord1[1]-coord2[1], 2) + pow(coord1[2]-coord2[2], 2))

#called by 'process'
#OBJ:  to check if it's covalently bound
def check_pro_lig_covalent(pro_chain, heterogens, lig_valid, conect_dict):
    #1. covalently bound by analyzing 'conect_dict'
    pro_atom_index = set([int(item[6:11]) for item in pro_chain])
    for i in xrange(len(lig_valid)):
        if lig_valid[i] != 1:
            continue
        covalent = False
        for atm in heterogens[i]:
            if not conect_dict.has_key(atm.atom_no):
                continue
            for index in conect_dict[atm.atom_no]:
                if index in pro_atom_index:
                    covalent = True
                    print "    ligand atom %d is covalently bound to protein atom %d"%(atm.atom_no, index)
                    break
            if covalent:
                break
        if covalent:
            lig_valid[i] = 0
    
    if 1 not in lig_valid:
        return
    #2. check steric clashes - maybe slow !!!!!!!!!!!!
    for line in pro_chain:
        if line.startswith("TER"):
            continue
        pro_coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        for i in xrange(len(lig_valid)):
            if lig_valid[i] != 1:
                continue
            clash = False
            for atm in heterogens[i]:
                dist = calc_distance(pro_coords, atm.coords)
                if dist <= clash_dist:
                    clash = True
                    print "    steric clash"
                    print "    distance between ligand atom %d and protein atom %s: %g"%\
                            (atm.atom_no, line[6:11].strip(), dist)
                    break
            if clash:
                lig_valid[i] = 0
        if 1 not in lig_valid:
            break


#called by 'process'
#OBJ: to check if any non-hydrogen atom on one ligand molecule 
#     was within 8.0 A of any non-hydrogen atom on another ligand molecule
#input:   heterogens    [<type Residue>, <type Polymer>, ...]
#         ligand_valid  [0,1,2,...]   (0 - invalid; 1 - valid; 2 - ion)
#ouput:   True(valid)/False(invalid)
def check_binary_complex(heterogens, ligand_valid):
    for i in xrange(len(ligand_valid)-1):
        if ligand_valid[i] != 1:
            continue
        for j in xrange(i+1, len(ligand_valid)):
            if ligand_valid[j] != 1:
                continue
            if heterogens[i].chain_id() != heterogens[j].chain_id():
                continue
            for atm_i in heterogens[i]:
                for atm_j in heterogens[j]:
                    dist = calc_distance(atm_i.coords, atm_j.coords)
                    if dist <= 8.0:
                        print "    ligands close to each other"
                        print "    distance between atom(%d) and atom(%d): %g"%(atm_i.atom_no, atm_j.atom_no, dist)
                        return False
    return True


#called by 'extract_pro'
def check_within(coord, heterogen):
    for atm in heterogen:
        dist = calc_distance(coord, atm.coords)
        if dist <= 14:
            return True
    return False

#called by 'extract_pro_ligand_ion_chain'
def extract_pro_chains(pro_chain, heterogen):
    #do nothing when the complex contains only one chain
    chain_ids = set([line[21] for line in pro_chain])
    if len(chain_ids) == 1:
        return [i for i in xrange(len(pro_chain))]
    #determine protein chains to be kept
    keep_chain_ids = set([])
    prev_chain = ""
    keep = False
    for i in xrange(len(pro_chain)):
        if pro_chain[i].startswith("TER"):
            continue
        chain = pro_chain[i][21]
        if chain != prev_chain:
            if keep and prev_chain!="":
                keep_chain_ids.add(prev_chain)
            prev_chain = chain
            keep = False
        if keep:
            continue
        x = float(pro_chain[i][30:38])
        y = float(pro_chain[i][38:46])
        z = float(pro_chain[i][46:54])
        if check_within((x,y,z), heterogen):
            keep = True
    if keep:
        keep_chain_ids.add(prev_chain)
    keep_pro_index = []
    keep = False
    for i in xrange(len(pro_chain)):
        if pro_chain[i].startswith("TER"):
            if keep:
                keep_pro_index.append(i)
            keep = False
            continue
        if pro_chain[i][21] in keep_chain_ids:
            keep_pro_index.append(i)
            keep = True
    return keep_pro_index


#OBJ:   to extract protein chain(s), unique ligand(s), and ion(s)
#       if there exist more than 1 unique ligand(s), a WARNING will be given!!
#input:  ligands     list of <type string>
#        heterogens  [<type Residue>, <type Polymer>, ...]
#        lig_valid   [0,1,2,...]   (0 - invalid; 1 - valid; 2 - ion)
#return: ...
def extract_pro_ligand_ion_chain(pro_chain, ligands, heterogens, lig_valid):
    keep_ions_list = []
    lig_name_chain = {}   #key=res_name, value=[(chain,index),...]
    for i in xrange(len(lig_valid)):
        if lig_valid[i] == 0:
            continue
        if isinstance(heterogens[i], Polymer):
            name = []
            chain = ""
            for res in heterogens[i].residues:
                name.append(res.res_name)
                chain = res.chain
            name.sort()
            name = "_".join(name)
            if not lig_name_chain.has_key(name):
                lig_name_chain[name] = []
            lig_name_chain[name].append((chain,i))
        elif isinstance(heterogens[i], Residue):
            if lig_valid[i] == 2:
                keep_ions_list.append("%s%4d"%(heterogens[i].chain, heterogens[i].res_no))
                continue
            if not lig_name_chain.has_key(heterogens[i].res_name):
                lig_name_chain[heterogens[i].res_name] = []
            lig_name_chain[heterogens[i].res_name].append((heterogens[i].chain,i))
        else:
            pass
    #only the first one of unique ligands were kept
    if len(lig_name_chain.keys()) > 1:
        print "    (MULTI-LIGANDS)contains %d (unique) valid ligands"%len(lig_name_chain.keys())
    keep_lig_list = []        #[[(chain,res_no),...], ...]   list of valid ligands
    keep_lig_index = []
    keep_pro_index = []       #[[line_no,...],...]  list of corresponding protein chain(s)
    keep_ions_index = []
    for key,value in lig_name_chain.iteritems():
        if len(value) > 1:
            print "    ligand(%s) appears %d times"%(key, len(value))
        i = value[0][1]
        each_keep_pro_index = extract_pro_chains(pro_chain, heterogens[i])
        keep_pro_index.append(each_keep_pro_index)
        keep_pro_chain_id = set([pro_chain[j][21] for j in each_keep_pro_index])
        keep_ions_index.append([j for j in xrange(len(ligands)) \
                if (ligands[j][21:26] in keep_ions_list and ligands[j][21] in keep_pro_chain_id)])
        each_keep_lig_list = []
        if isinstance(heterogens[i], Polymer):
            for res in heterogens[i].residues:
                each_keep_lig_list.append("%s%4d"%(res.chain, res.res_no))
        else:
            each_keep_lig_list.append("%s%4d"%(heterogens[i].chain, heterogens[i].res_no))
        each_keep_lig_index = [j for j in xrange(len(ligands)) if ligands[j][21:26] in each_keep_lig_list]
        keep_lig_index.append(each_keep_lig_index)
    
    return keep_pro_index,keep_lig_index,keep_ions_index


#called by 'process'
#OBJ: to output protein chains, ions and waters
def output_protein(outf_pro, pro_chain, keep_pro_index, ligands, keep_ions_index, conect):
    atm_index = set([])
    for i in keep_pro_index:
        outf_pro.write(pro_chain[i])
        if not pro_chain[i].startswith("TER"):
            atm_index.add(int(pro_chain[i][6:11]))
    for i in keep_ions_index:
        outf_pro.write(ligands[i])
        atm_index.add(int(ligands[i][6:11]))
    for line in conect:
        i1 = line[6:11]; i2 = line[11:16]
        i3 = line[16:21].strip(); i4 = line[21:26].strip()
        i5 = line[26:31].strip()
        if int(i1) not in atm_index:
            continue
        if int(i2) not in atm_index:
            continue
        if i3!="" and (int(i3) not in atm_index):
            continue
        if i4!="" and (int(i4) not in atm_index):
            continue
        if i5!="" and (int(i5) not in atm_index):
            continue
        outf_pro.write(line)

#called by 'process'
#OBJ:  to output valid ligands
def output_ligands(outf_lig, ligands, keep_lig_index, conect):
    atm_index = set([])
    for i in keep_lig_index:
        outf_lig.write(ligands[i])
        atm_index.add(int(ligands[i][6:11]))
    for line in conect:
        i1 = line[6:11]; i2 = line[11:16]
        i3 = line[16:21].strip(); i4 = line[21:26].strip()
        i5 = line[26:31].strip()
        if int(i1) not in atm_index:
            continue
        if int(i2) not in atm_index:
            continue
        if i3!="" and (int(i3) not in atm_index):
            continue
        if i4!="" and (int(i4) not in atm_index):
            continue
        if i5!="" and (int(i5) not in atm_index):
            continue
        outf_lig.write(line)

def output_waters(outf_water, waters, keep_water_index, conect):
    atm_index = set([])
    for i in keep_water_index:
        outf_water.write(waters[i])
        atm_index.add(int(waters[i][6:11]))
    for line in conect:
        i1 = line[6:11]; i2 = line[11:16]
        i3 = line[16:21].strip(); i4 = line[21:26].strip()
        i5 = line[26:31].strip()
        if int(i1) not in atm_index:
            continue
        if int(i2) not in atm_index:
            continue
        if i3!="" and (int(i3) not in atm_index):
            continue
        if i4!="" and (int(i4) not in atm_index):
            continue
        if i5!="" and (int(i5) not in atm_index):
            continue
        outf_water.write(line)


def process(pdb_file, out_dir, outf_log=None):
    print "\n\n>>> to process",pdb_file
    inf = open(pdb_file,"r")
    #1. check resolution
    line = inf.readline()
    while line!="" and not line.startswith("REMARK   2 RESOLUTION."):
        line = inf.readline()
    if line != "":
        #resolution = float(line[23:30])
        temp = line.split()[-2]
        resolution = float(temp)
        #print "  =>resolution=%s"%(line[23:30].strip())
        print "  =>resolution=%s"%(line.split()[-2])
        if resolution > 3.0:
            print "  =>(fail) resolution",line[23:30].strip(),"is larger than 3.0 A"
            return False
    #2. read PDB file and check if protein contains any DNA or RNA
    inf.seek(0)
    pro_chain,ligands,waters,conect,conect_dict = read_pdb_file(inf)
    inf.close()
    if len(pro_chain) == 0:
        print "  =>(fail) can't find protein chain(s)"
        return False
    if len(ligands) == 0:
        print "  =>(fail) can't find ligand(s)"
        return False
    #3. construct 'Residue' or 'Polymer'
    print "  STEP 1. to call construct_ligands"
    heterogens = construct_ligands(ligands, conect_dict)
    print "    totally %d candidate heterogens (ligand, solvent, metal, ion etc.)"%len(heterogens)
    for i in xrange(len(heterogens)):
        print "    %-2d: %s"%(i+1, heterogens[i].get_name())
    #4. check if contains valid ligand
    print "  STEP 2. to call check_valid_ligands"
    lig_valid = check_valid_ligand(heterogens)
    if 1 not in lig_valid:
        #print "  =>(fail) contains no valid ligand(s) or special ligand(s)"
        print "  =>(fail) no more valid ligand"
        return False
    #4.2 check if ligand is covalently bound
    print "  STEP 3. to call check_pro_lig_covalent/steric clash"
    check_pro_lig_covalent(pro_chain, heterogens, lig_valid, conect_dict)
    if 1 not in lig_valid:
        #print "  =>(fail) contain covalently bound ligand(s) or steric clash"
        print "  =>(fail) no more valid ligand"
        return False
    #5.check if it's binary complex
    print "  STEP 4. to call check_binary_complex"
    is_valid = check_binary_complex(heterogens, lig_valid)
    if not is_valid:
        #print "  =>(fail) contain ligands that close to each other"
        print "  =>(fail) no more valid ligand"
        return False
    """
    #6. check protein (contains any non-standard residues)
    error = check_protein(pro_chain, ligands, lig_valid)
    if error:
        #more detailed information, which residue...
        print pdb_file+": contains non-standard residue"
        return None
    #7. check ligands with low buried surface area (<= 15%)
    """
    #6. extract protein chains, valid ligand(s) and ions if any.
    #   for protein chains, only those within 14 A of valid ligand(s) were kept
    print "  STEP 5. to call extract_pro_ligand_ion_chain"
    keep_pro_index,keep_lig_index,keep_ions_index = \
            extract_pro_ligand_ion_chain(pro_chain, ligands, heterogens, lig_valid)
    #7. output proteins, valid ligand(s), ions and waters
    print "  STEP 6. to save protein and valid ligand(s)"
    print "    totally %d valid ligand(s) being kept!"%len(keep_lig_index)
    base_name = os.path.basename(pdb_file)
    base_name = base_name[:base_name.rfind(".")]
    if len(keep_pro_index) == 1:
        pro_file = os.path.join(out_dir, base_name+"_pro.pdb")
        lig_file = os.path.join(out_dir, base_name+"_lig.pdb")
        if outf_log is not None:
            print >>outf_log, pro_file,lig_file
        #water_file = os.path.join(out_dir, base_name+"_water.pdb")
        outf_pro = open(pro_file,"w")
        outf_lig = open(lig_file, "w")
        #outf_water = open(water_file, "w")
        #output protein and ligands
        output_protein(outf_pro, pro_chain, keep_pro_index[0], ligands, keep_ions_index[0], conect)
        print >>outf_pro, "END"
        outf_pro.close()
        output_ligands(outf_lig, ligands, keep_lig_index[0], conect)
        print >>outf_lig, "END"
        outf_lig.close()
    else:
        for i in xrange(len(keep_pro_index)):
            pro_file = os.path.join(out_dir, base_name+"_pro_%d.pdb"%(i+1))
            lig_file = os.path.join(out_dir, base_name+"_lig_%d.pdb"%(i+1))
            if outf_log is not None:
                print >>outf_log, pro_file, lig_file
            outf_pro = open(pro_file, "w")
            outf_lig = open(lig_file, "w")
            output_protein(outf_pro, pro_chain, keep_pro_index[i], ligands, keep_ions_index[i], conect)
            print >>outf_pro, "END"
            outf_pro.close()
            output_ligands(outf_lig, ligands, keep_lig_index[i], conect)
            print >>outf_lig, "END"
            outf_lig.close()
    #output water - only those having same chain ID with kept protein
    """
    keep_pro_chain_id = set([pro_chain[i][21] for i in keep_pro_index])
    keep_water_index = [i for i in xrange(len(waters)) if waters[i][21] in keep_pro_chain_id]
    output_waters(outf_water, waters, keep_water_index, conect)
    print >>outf_water, "END"
    outf_water.close()
    """
    print "  =>(success)"

    return True

