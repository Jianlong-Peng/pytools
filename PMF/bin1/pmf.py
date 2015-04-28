'''
#=============================================================================
#     FileName: pmf.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-03-22 11:06:21
#   LastChange: 2014-07-07 11:02:51
#      History:
#=============================================================================
'''
import sys
import math
import pybel
import atom_type

Rmax = 12  #A
dr = 0.2   #A
min_num_atom_pair = 1000
KB = 0.00198721  #Kcal/mol/k
T = 298  #k

def dist2idx_v1(dist):
    '''
    if dist < 4, return 0
    else, return int(math.floor(dist - 3))
    '''
    assert dist < Rmax
    index = 0
    if dist >= 4.:
        index = int(math.floor(dist - 3))
    return index

dist2idx_v2 = lambda dist: int(math.ceil(dist/dr))-1

def get_pair_index(_type, ptype, ltype, dist=None):
    '''
    _type: which kind of atom type being used
      0 - PMF04 atom type
    ptype: protein atom type
    ltype: ligand atom type
    dist : None or float indicates distance between protein and ligand atom
    '''
    idx_row = pair2index(_type, ptype, ltype)
    if dist is None:
        idx_col = None
    else:
        idx_col = dist2idx_v2(dist)
    return (idx_row, idx_col)

def index2pair(_type, index):
    '''
    _type: which kind of atom type being used
      0 - PMF04 atom type
    index: int, [0, PTYPE*LTYPE]
    '''
    pidx = index / len(atom_type.LTYPE[_type])
    lidx = index % len(atom_type.LTYPE[_type])
    return (atom_type.PTYPE[_type][pidx], atom_type.LTYPE[_type][lidx])

def pair2index(_type, ptype, ltype):
    '''
    _type: which kind of atom type being used
      0 - PMF04 atom type
    ptype: protein atom type
    ltype: ligand atom type
    '''
    pidx = atom_type.PTYPE[_type].index(ptype)
    lidx = atom_type.LTYPE[_type].index(ltype)
    index = pidx * len(atom_type.LTYPE[_type])+ lidx
    return index

class NumAtomPairs:
    def __init__(self):
        '''
        nij: number of atom pairs ij (protein vs ligand) in spherical shell
             (len(PTYPE)*len(LTYPE)) X int(Rmax/dr)
        nkj: number of pairs of protein atoms k of any type around a ligand atom of type j
             (len(LTYPE)) X 9 (Rmax-3)
        nlj: number of pairs of ligand atoms l of any type around a ligand atom of type j
             (len(LTYPE)) X 9 (Rmax-3)
        (summation over all complexes)
        _type: int, atom type being used
          0 - PMF04 atom type

        nij is generated in each spherical shells, and the bin size is set to be 0.2 A.
        nkj and nlj are generated in spherical shells [0,4), [4,5), [5,6), ..., [11,12)
        they will be used to calculate volume correction factor
        '''
        self._type = -1
        #number of atom pairs ij. (PTYPE*LTYPE) rows X (Rmax/dr) cols
        self.nij = [] #[[0 for i in xrange(num_bins)] for j in xrange(num_ptypes*num_ltypes)]
        #number of pairs of protein atoms k of any type around a ligand atom of type j
        #(LTYPE) rows X 9 cols
        self.nkj = [] #[[0 for i in xrange(int(Rmax-3))] for j in xrange(num_ltypes)]
        #number of pairs of ligand atoms l of any type around a ligand atom of type j
        #(LTYPE) rows X 9 cols
        self.nlj = [] #[[0 for i in xrange(int(Rmax-3))] for j in xrange(num_ltypes)]

    def __str__(self):
        '''
        to display N_ij (protein ligand atom pair ij)
        '''
        lines = 'atom type %d\n'%self._type
        for i in xrange(len(self.nij)):
            ptype,ltype = index2pair(self._type, i)
            lines += 'N_ij (pro_lig) %s_%s\n'%(ptype,ltype)
            for j in xrange(len(self.nij[i])):
                lines += ' %7.1f'%((j+1)*dr)
            lines += '\n'
            for j in xrange(len(self.nij[i])):
                lines += ' %7d'%self.nij[i][j]
            lines += '\n'
            lines += 'N_ij_bulk %d\n'%sum(self.nij[i])
            lines += '\n'
        return lines

    def get_nij(self, ptype, ltype, dist=None):
        '''
        to get 'nij' according to atom type 'ptype' and 'ltype'.
        if dist is float, distance will be mapped to bin, and return 'nij(r)',
        otherwise, return list (number of atom pairs ij in all bins)
        '''
        idx_row, idx_col = get_pair_index(self._type, ptype, ltype, dist)
        if idx_col is None:
            return self.nij[idx_row]
        else:
            return self.nij[idx_row][idx_col]

    def get_nkj(self, ltype, dist=None):
        '''
        to get 'nkj' according to atom type 'ltype'.
        if 'dist' is float, distance will be mapped to bin, and return 'nkj(r)',
        otherwise, return list (number of atom pairs jk in all the bins)
        '''
        idx_row = atom_type.LTYPE[self._type].index(ltype)
        if dist is None:
            return self.nkj[idx_row]
        else:
            idx_col = dist2idx_v1(dist)
            return self.nkj[idx_row][idx_col]

    def get_nlj(self, ltype, dist=None):
        '''
        to get 'nlj' according to atom type 'ltype'.
        if 'dist' is float, distance will be mapped to bin, and return 'nlj(r)',
        otherwise, return list (number of atom pairs lj in all bins)
        '''
        idx_row = atom_type.LTYPE[self._type].index(ltype)
        if dist is None:
            return self.nlj[idx_row]
        else:
            idx_col = dist2idx_v1(dist)
            return self.nlj[idx_row][idx_col]

    def update(self, infile):
        '''
        to extract atom pair information from lines starts with 'PL' or 'LL',
        and update 'self.nij', 'self.nkj', and 'self.nlj'.
        'infile' should be the one generated by 'pmf_step1.py'
        '''
        count_pl = 0
        count_ll = 0
        inf = open(infile,"r")
        line = inf.readline()
        assert line.startswith("atom type")
        tmp_type = int(line.split()[-1])
        if self._type == -1:
            self._type = tmp_type
            num_bins = int(Rmax/dr)
            num_ptypes = len(atom_type.PTYPE[self._type])
            num_ltypes = len(atom_type.LTYPE[self._type])
            self.nij = [[0 for i in xrange(num_bins)] for j in xrange(num_ptypes*num_ltypes)]
            self.nkj = [[0 for i in xrange(int(Rmax-3))] for j in xrange(num_ltypes)]
            self.nlj = [[0 for i in xrange(int(Rmax-3))] for j in xrange(num_ltypes)]
        else:
            if tmp_type != self._type:
                print "Error: incompatible atom type found in %s"%infile
                print "       already(%d) VS. new(%d)"%(self._type, tmp_type)
                sys.exit(1)
        for line in inf:
            if line.startswith("PL"):
                temp = line.strip().split(",")
                #update self.nij
                dist = float(temp[-1])
                if dist >= Rmax:
                    print "Warning: distance is larger than %g"%Rmax
                    sys.stdout.write(line)
                    continue
                idx_row, idx_col = get_pair_index(self._type, temp[1], temp[2], dist)
                self.nij[idx_row][idx_col] += 1
                #update self.nkj
                idx_row = atom_type.LTYPE[self._type].index(temp[2])
                idx_col = dist2idx_v1(dist)
                self.nkj[idx_row][idx_col] += 1
                count_pl += 1
            elif line.startswith("LL"):
                temp = line.strip().split(",")
                #update self.nlj
                idx_row = atom_type.LTYPE[self._type].index(temp[1])
                dist = float(temp[-1])
                if dist >= Rmax:
                    print "Warning: distance is larger than %g"%Rmax
                    sys.stdout.write(line)
                    continue
                idx_col = dist2idx_v1(dist)
                self.nlj[idx_row][idx_col] += 1
                count_ll += 1
            else:
                pass
        inf.close()
        print "totally read %d protein ligand atom pairs and %d ligand ligand atom pairs from file %s"%(count_pl,count_ll,infile)

    def write(self, outfile):
        '''
        to write 'nij', 'nkj', and 'nlj' to file 'outfile'
        '''
        if self._type == -1:
            print "Warning: <NumAtomPairs> not being initialized!"
            sys.exit(1)
        outf = open(outfile,"w")
        outf.write("atom type %d\n"%self._type)
        outf.write("Nij (pro_lig)\n")
        for i in xrange(len(self.nij)):
            ptype,ltype = index2pair(self._type, i)
            outf.write(ptype+"_"+ltype)
            for j in xrange(len(self.nij[i])):
                outf.write(" %d"%self.nij[i][j])
            outf.write("\n")
        outf.write("\nNkj (lig)\n")
        for i in xrange(len(self.nkj)):
            ltype = atom_type.LTYPE[self._type][i]
            outf.write(ltype)
            for j in xrange(len(self.nkj[i])):
                outf.write(" %d"%self.nkj[i][j])
            outf.write("\n")
        outf.write("\nNlj (lig)\n")
        for i in xrange(len(self.nlj)):
            ltype = atom_type.LTYPE[self._type][i]
            outf.write(ltype)
            for j in xrange(len(self.nlj[i])):
                outf.write(" %d"%self.nlj[i][j])
            outf.write("\n")
        outf.write("\n")
        outf.close()

    def load(self, infile):
        '''
        'infile' should be the one generated by 'self.write'
        'nij','nkj',and 'nlj' will be replaced!!!!!!!
        '''
        inf = open(infile,"r")
        #read self._type
        line = inf.readline()
        assert line.startswith("atom type")
        self._type = int(line.split()[-1])
        num_bins = int(Rmax/dr)
        num_ptypes = len(atom_type.PTYPE[self._type])
        num_ltypes = len(atom_type.LTYPE[self._type])
        self.nij = [[0 for i in xrange(num_bins)] for j in xrange(num_ptypes*num_ltypes)]
        self.nkj = [[0 for i in xrange(int(Rmax-3))] for j in xrange(num_ltypes)]
        self.nlj = [[0 for i in xrange(int(Rmax-3))] for j in xrange(num_ltypes)]
        #read Nij
        line = inf.readline()
        assert line.strip() == "Nij (pro_lig)"
        line = inf.readline()
        while line!="" and line.strip()!="":
            line = line.split()
            ptype,ltype = line[0].split("_")
            index = pair2index(self._type, ptype, ltype)
            self.nij[index] = map(int, line[1:])
            line = inf.readline()
        #read Nkj
        line = inf.readline()
        assert line.strip() == "Nkj (lig)"
        line = inf.readline()
        while line!="" and line.strip()!="":
            line = line.split()
            index = atom_type.LTYPE[self._type].index(line[0])
            self.nkj[index] = map(int, line[1:])
            line = inf.readline()
        #read Nlj
        line = inf.readline()
        assert line.strip() == "Nlj (lig)"
        line = inf.readline()
        while line!="" and line.strip()!="":
            line = line.split()
            index = atom_type.LTYPE[self._type].index(line[0])
            self.nlj[index] = map(int, line[1:])
            line = inf.readline()
        #done
        inf.close()


def smooth(fj):
    '''
    - strategy:
    define set of spherical shells called segments seg(r)
    (seg = 1,2,...,R/m) with thickness m consecutively
    separating our sphere with radius R into R/m segments.
    Here, m was chosen to be 0.2 A
    the smoothed volume correction factor will be:
    f_{j}'(seg) = \frac{1}{17}*\sum_{i=seg-8}^{seg+8} f_{j}(seg)
    '''
    u = []
    d = 0.
    while d < Rmax:
        j = dist2idx_v1(d)
        u.append(fj[j])
        d += dr
    w = []
    for j in xrange(len(u)):
        val = u[j]
        if j>=8 and j<len(u)-8:
            val = 0.
            for k in xrange(17):
                val += u[j-8+k]
            val /= 17.
        w.append(val)
    return w


class Potential:
    def __init__(self):
        '''
        self._type: atom type being used
          0 - PMF04 atom type
        self.Aij: without volume correction factor fj(r)
        self.Aij_f: use volume correction factor fj(r)
        both (PTYPE*LEYPT) rows X (Rmax/dr) cols
        '''
        self._type = -1
        self.Aij = []
        self.Aij_f = []

    def construct(self, ap):
        '''
        to calculate pair potential A_{ij}(r) from numbers of atom pairs ij.

        A_{ij}(r) = -K_{B}*T*ln{f_{j}(r)*\frac{rho_{ij}(r)}{rho_{ij,bulk}}}
        where, f_{j}(r) is the volume correction factor
               rho_{ij}(r) = \sum_{m}\frac{n_{ij,m}(r)}{4*pi*r^2*dr}
               rho_{ij,bulk} = \sum_{m}\frac{N_{ij,m}}{4*pi*R^3/3}
        
        ap: <type NumAtomPairs>
        '''
        self._type = ap._type
        num_ptypes = len(atom_type.PTYPE[self._type])
        num_ltypes = len(atom_type.LTYPE[self._type])
        #1. calculate volume correction factor
        fj = []
        for i in xrange(num_ltypes):
            nkj_bulk = sum(ap.nkj[i])
            nlj_bulk = sum(ap.nlj[i])
            vj_bulk = 0.
            if (nkj_bulk+nlj_bulk) > 0:
                vj_bulk = float(nkj_bulk) / (nkj_bulk + nlj_bulk)
            nbins = int(Rmax-3)
            each_fj = [1.0 for j in xrange(nbins)]
            for j in xrange(nbins):
                nkj = ap.nkj[i][j]
                nlj = ap.nlj[i][j]
                vj = 0.
                if (nkj+nlj) > 0:
                    vj = float(nkj)/(nkj+nlj)
                #if vj_bulk > 1E-5:
                #    each_fj[j] = vj / vj_bulk
                if vj > 1E-5:
                    each_fj[j] = vj_bulk / vj
            fj.append(smooth(each_fj))
        #2. calculate protein-ligand atom pair potential A_ij(r)
        vol_bulk = 4 * math.pi / 3. * (Rmax**3)
        for i in xrange(num_ptypes*num_ltypes):
            lidx = i % num_ltypes
            nij_bulk = sum(ap.nij[i])
            w = []
            w_f = []
            for j in xrange(len(ap.nij[i])):
                nij = ap.nij[i][j]
                r1 = dr * j
                r2 = dr * (j+1)
                vol = 4 * math.pi / 3. * (r2**3 - r1**3)
                rho_ij = nij / vol
                rho_bulk = nij_bulk / vol_bulk

                if nij_bulk < min_num_atom_pair:
                    each_Aij = 0.
                    each_Aij_f = 0.
                elif float(nij)/nij_bulk < 1E-5:
                    each_Aij = 3.
                    each_Aij_f = 3.
                else:
                    each_Aij = KB*T*(math.log(rho_bulk) - math.log(rho_ij))
                    each_Aij_f = KB*T*(math.log(rho_bulk) - math.log(rho_ij) - math.log(fj[lidx][j]))
                w.append(each_Aij)
                w_f.append(each_Aij_f)
            self.Aij.append(w)
            self.Aij_f.append(w_f)

    def __str__(self):
        lines = 'atom type %d\n'%self._type
        for i in xrange(len(self.Aij)):
            ptype,ltype = index2pair(self._type, i)
            lines += 'atom pair (pro_lig) %s_%s\n'%(ptype, ltype)
            for j in xrange(len(self.Aij[i])):
                lines += ' %7.1f'%(dr*j)
            lines += '\n'
            for Aij in self.Aij[i]:
                lines += ' %7.6f'%Aij
            lines += '\n'
            for Aij_f in self.Aij_f[i]:
                lines += ' %7.6f'%Aij
            lines += '\n\n'
        return lines
    
    def get_Aij(self,ptype, ltype, dist=None):
        idx_row, idx_col = get_pair_index(self._type, ptype, ltype, dist)
        if idx_col is None:
            return self.Aij[idx_row]
        else:
            return self.Aij[idx_row][idx_col]
    
    def get_Aij_f(self,ptype, ltype, dist=None):
        idx_row, idx_col = get_pair_index(self._type, ptype, ltype, dist)
        if idx_col is None:
            return self.Aij_f[idx_row]
        else:
            return self.Aij_f[idx_row][idx_col]

    def write(self, outfile):
        outf = open(outfile,"w")
        outf.write("atom type %d\n"%self._type)
        outf.write("Aij (pro_lig)\n")
        for i in xrange(len(self.Aij)):
            ptype,ltype = index2pair(self._type, i)
            outf.write(ptype+"_"+ltype)
            for j in xrange(len(self.Aij[i])):
                outf.write(" %.8f"%self.Aij[i][j])
            outf.write("\n")
        outf.write("\nAij_f (pro_lig)\n")
        for i in xrange(len(self.Aij_f)):
            ptype,ltype = index2pair(self._type, i)
            outf.write(ptype+"_"+ltype)
            for j in xrange(len(self.Aij_f[i])):
                outf.write(" %.8f"%self.Aij_f[i][j])
            outf.write("\n")
        outf.write("\n")
        outf.close()

    def load(self, infile):
        '''
        infile: should be the one generated by 'self.write'
        '''
        assert len(self.Aij) == 0
        inf = open(infile,"r")
        #read self._type
        line = inf.readline()
        assert line.startswith("atom type")
        self._type = int(line.split()[-1])
        num_ptypes = len(atom_type.PTYPE[self._type])
        num_ltypes = len(atom_type.LTYPE[self._type])
        #read Aij
        line = inf.readline()
        assert line.strip() == "Aij (pro_lig)"
        self.Aij = [None for i in xrange(num_ptypes*num_ltypes)]
        line = inf.readline()
        while line!="" and line.strip()!="":
            line = line.split()
            ptype,ltype = line[0].split("_")
            i = pair2index(self._type, ptype, ltype)
            self.Aij[i] = map(float, line[1:])
            line = inf.readline()
        #read Aij_f
        line = inf.readline()
        assert line.strip() == "Aij_f (pro_lig)"
        self.Aij_f = [None for i in xrange(num_ptypes*num_ltypes)]
        line = inf.readline()
        while line!="" and line.strip()!="":
            line = line.split()
            ptype,ltype = line[0].split("_")
            i = pair2index(self._type, ptype, ltype)
            self.Aij_f[i] = map(float, line[1:])
            line = inf.readline()
        inf.close()


#OBJ: to cunt the number of atom pair ij in the spherical shell and the reference sphere
#     for the complex (protein, ligand)
#INPUT : 
#      protein: OBMol
#      ligand : OBMol
#      ptypes : None, or list of atom types of protein generated by atom_type.parse_protein
#RETURN: dict, key="ptype_ltype", value=[occurance,...]
def count_number_of_pairs(protein, ligand, ptypes=None):
    assert False
    ltypes = atom_type.parse_ligand(ligand)
    pmf_number = {}
    for p in atom_type.PTYPE[atom_type._type]:
        for l in atom_type.LTYPE[atom_type._type]:
            pmf_number[p+"_"+l] = [0 for _ in xrange(int(Rmax/dr))]
    pidx = -1
    for patom in pybel.ob.OBMolAtomIter(protein):
        pidx += 1
        if ptypes is None:
            ptype = atom_type.parse_protein_atom[atom_type._type](patom)
        else:
            ptype = ptypes[pidx]
        if ptype == "":
            continue
        lidx = -1
        for latom in pybel.ob.OBMolAtomIter(ligand):
            lidx += 1
            ltype = ltypes[lidx]
            if ltype == "":
                continue
            dist = patom.GetDistance(latom)
            if dist > Rmax:
                continue
            key = ptype + "_" + ltypes[lidx]
            index = int(math.ceil(dist/dr)) - 1
            pmf_number[key][index] += 1    #key must be in pmf_number

    return pmf_number

#OBJ: to parse atoms of protein and ligand, and return the atom types and distance
#INPUT:
#  protein: OBMol
#  ligand : OBMol
#  _type  : specify the atom type set to be used
#  ptypes : None or list of atom types of protein generated by atom_type.parse_protein
#  ltypes : None or list of atom types of ligand generated by atom_type.parse_ligand
#RETURN
#  ["PL,ptype,ltype,pidx,lidx,dist\n",...]
def gen_atom_pair_PL(protein, ligand, _type, ptypes=None, ltypes=None):
    if ltypes is None:
        ltypes = atom_type.parse_ligand(ligand, _type)
    if ptypes is None:
        ptypes = atom_type.parse_protein(protein, _type)
    result = []
    for patom in pybel.ob.OBMolAtomIter(protein):
        pidx = patom.GetIdx()-1
        ptype = ptypes[pidx]
        if ptype == "":
            continue
        for latom in pybel.ob.OBMolAtomIter(ligand):
            lidx = latom.GetIdx()-1
            ltype = ltypes[lidx]
            if ltype == "":
                continue
            dist = patom.GetDistance(latom)
            if dist >= Rmax:
                continue
            each_result = "PL,%s,%s,%d,%d,%.6f\n"%(ptype,ltype,pidx+1,lidx+1,dist)
            result.append(each_result)
    return result

#OBJ: return atom pairs of ligand atom i and j, as well as distance
#INPUT:
#  ligand: OBMol
#  _type : specify the atom type set to be used
#  ltypes: None or list of atom types of protein generated by atom_type.parse_ligand
#RETURN:
#  ["LL,ltype,ltype,lidx,lidx,dist\n",...]
def gen_atom_pair_LL(ligand, _type, ltypes=None):
    if ltypes is None:
        ltypes = atom_type.parse_ligand(ligand, _type)
    result = []
    for latom1 in pybel.ob.OBMolAtomIter(ligand):
        lidx1 = latom1.GetIdx()-1
        if ltypes[lidx1] == "":
            continue
        for latom2 in pybel.ob.OBMolAtomIter(ligand):
            lidx2 = latom2.GetIdx()-1
            if lidx2==lidx1 or ltypes[lidx2]=="":
                continue
            dist = latom1.GetDistance(latom2)
            if dist >= Rmax:
                continue
            each_result = "LL,%s,%s,%d,%d,%.6f\n"%(ltypes[lidx1],ltypes[lidx2],lidx1+1,lidx2+1,dist)
            result.append(each_result)
    return result

#INPUT: (dict, FILE)
def write_number_of_pairs(pmf_number, outf):
    keys = pmf_number.keys()
    keys.sort()
    for key in keys:
        outf.write(key)
        for value in pmf_number[key]:
            outf.write(" %d"%value)
        outf.write("\n")


