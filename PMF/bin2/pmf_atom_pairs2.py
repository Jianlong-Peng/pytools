'''
#=============================================================================
#     FileName: pmf_atom_pairs2.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-08-14 19:14:37
#   LastChange: 2014-08-14 13:14:43
#      History:
#=============================================================================
'''
import sys
from glob import glob
import threading
import Queue
import pmf_atom_typer
import pybel

global mutex, ligand, protein, outf, ap_iter

class AtomPairWalker(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)

    def run(self):
        global mutex, ligand, protein, outf, ap_iter
        while True:
            mutex.acquire()
            try:
                i,j = ap_iter.next()
            except StopIteration:
                break
            mutex.release()
            dist = protein.atoms[j].OBAtom.GetDistance(ligand.atoms[i].OBAtom)
            if dist >= 12.0:
                continue 
            ltype = pmf_atom_typer.LigAtomTyper(ligand.atoms[i].OBAtom)
            if ltype == "":
                continue
            ptype = pmf_atom_typer.ProAtomTyper(protein.atoms[j].OBAtom)
            if ptype in ("","HD"):
                continue
            mutex.acquire()
            print >>outf, "%s,%d,%s,%d,%.6f"%(ltype, i+1, ptype, j+1, dist)
            #outf.flush()
            mutex.release()


class AtomPairs:
    def __init__(self):
        global ligand,protein
        self.max_latoms = len(ligand.atoms)
        self.max_patoms = len(protein.atoms)
        self.i = 0
        self.j = -1
        self.incre_i = False
    
    def next(self):
        global ligand,protein
        self.j += 1
        if self.j == self.max_patoms:
            self.j = 0
            self.incre_i = True
        while self.j<self.max_patoms and protein.atoms[self.j].OBAtom.IsHydrogen():
            self.j += 1
        if self.j == self.max_patoms:
            raise StopIteration
        if self.incre_i:
            self.i += 1
            self.incre_i = False
            while self.i<self.max_latoms and ligand.atoms[self.i].OBAtom.IsHydrogen():
                self.i += 1
            if self.i == self.max_latoms:
                raise StopIteration
        return (self.i, self.j)

    def __iter__(self):
        return self


def make_pairs():
    global ligand,protein,pair_queue
    pair_queue = Queue.Queue(0)
    for i in xrange(len(ligand.atoms)):
        if ligand.atoms[i].OBAtom.IsHydrogen():
            continue
        for j in xrange(len(protein.atoms)):
            if protein.atoms[i].OBAtom.IsHydrogen():
                continue
            if protein.atoms[i].OBAtom.GetResidue().GetName() in ('HOH','WAT'):
                continue
            pair_queue.put((i,j))

def do_for_each(name, num_threads):
    _format = name[name.rfind(".")+1:]
    basename = name[:name.rfind(".")]
    out_name = basename + ".num"
    pdb_id = basename.split("_")[-1]
    pdb_name = "/home/xmluo/jlpeng/cMet/pdb_protein/%s_protein.pdb"%pdb_id
    global ligand,protein,outf,ap_iter
    ligand = pybel.readfile(_format, name).next()
    protein = pybel.readfile("pdb", pdb_name).next()
    outf = open(out_name, "w")
    #make_pairs()
    #print "after creating pair_queue"
    ap_iter = AtomPairs()
    threads = []
    for i in xrange(num_threads):
        threads.append(AtomPairWalker())
    print "len(threads)=%d"%len(threads)
    for t in threads:
        t.start()
    for t in threads:
        t.join()
        if t.isAlive():
            print "thread %s is alive"%t.getName()
        else:
            print "thread %s is dead"%t.getName()
    outf.close()
    


def main(argv=sys.argv):
    if len(argv) < 2:
        print "\n  Usage: %s mol_in[...]"%argv[0]
        print "  mol_in: should be like cmet_xxx_{pdbid}.pdb"
        print "\n  information of atom pair for each molecule"
        print "  will be saved in a file named basename{mol_in}.num"
        print ""
        sys.exit(1)

    global mutex
    mutex = threading.Lock()

    for name in argv[1:]:
        for each_name in glob(name):
            do_for_each(each_name, 10)
main()

