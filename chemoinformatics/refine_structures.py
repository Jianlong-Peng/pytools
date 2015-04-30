'''
#=============================================================================
#     FileName: refine_structures.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2015-04-30 14:10:55
#   LastChange: 2015-04-30 16:47:29
#      History:
#=============================================================================
'''
import sys
from getopt import getopt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.SaltRemover import SaltRemover


def main(argv=sys.argv):
    valid_elements = ['H','C','N','O','F','Si','P','S','Cl','Br','I']
    valid_atomic_num = [1,6,7,8,9,14,15,16,17,35,53]

    if len(argv) < 3:
        print """
OBJ
  to strip self-defined counterions

Usage:
  %s [options] input output

[options]
  --strip            : if given, to run stripping salts/solvents
  --strip-sdf    file: specify the mols/fragments to be removed
  --strip-smarts file: one SMARTS string per line
  --filter-invalid   : if given, to remove molecules containing
                       R group or elements other than
                       %s
  --addh             : if given, to add hydrogens
  --make3d           : if given, 3D coordinates will be generated

Attention
  1. rdkit.Chem.SaltRemover.SaltRemover is called
  2. if neither `--strip-sdf` nor `--strip-smarts` is provided,
     stripping salts will be done according to default salts defined
     in `RDConfig.RDDataDir/Salts.txt`
  3. both `input` and `output` are .sdf
  4. whenever `--make3d` is given, please make sure that
     there is no complex or salts/solvents can be stripped.
     Otherwise, maybe there is something wrong with optimized structure.
"""%(argv[0], str(valid_elements))
        sys.exit(1)

    options,args = getopt(argv[1:],'',['strip-sdf=','strip-smarts=','strip','filter-invalid','addh','make3d'])
    filter_invalid = False
    strip = False
    addh = False
    make3d = False
    strip_sdf = None
    strip_smarts = None
    for opt,val in options:
        if opt == '--strip':
            strip = True
        elif opt == "--strip-sdf":
            strip_sdf = val
        elif opt == '--strip-smarts':
            strip_smarts = val
        elif opt == '--filter-invalid':
            filter_invalid = True
        elif opt == '--addh':
            addh = True
        elif opt == '--make3d':
            make3d = True
        else:
            print "Error: invalid option",opt
            sys.exit(1)
    assert len(args) == 2
    infile = args[0]
    outfile = args[1]

    smarts = ""
    if strip_sdf is not None:
        print "To load fragments from",strip_sdf
        count = 0
        for m in Chem.SDMolSupplier(strip_sdf):
            count += 1
            if m is None:
                print "Warning: failed to read %dth molecule in %s"%(count, strip_sdf)
                continue
            smarts += (Chem.MolToSmarts(m) + "\n")
    if strip_smarts is not None:
        print "to load fragments from",strip_smarts
        for line in open(strip_smarts,'r'):
            smarts += line

    if strip:
        if smarts == "":
            remover = SaltRemover(defnData=smarts)
        else:
            remover = SaltRemover()
    else:
        remover = None

    inf = Chem.SDMolSupplier(infile)
    outf = Chem.SDWriter(outfile)
    count = 0
    for m in inf:
        count += 1
        if m is None:
            print "Warning: failed to load %dth molecule from %s"%(count, infile)
            continue
        if filter_invalid:
            invalid = False
            for a in m.GetAtoms():
                if a.GetAtomicNum() not in valid_atomic_num:
                    invalid = True
                    break
            if invalid:
                continue
        if strip:
            m = remover.StripMol(m)
        if addh:
            m = AllChem.AddHs(m)
        if make3d:
            AllChem.EmbedMolecule(m)
            AllChem.UFFOptimizeMolecule(m)
        outf.write(m)
    outf.close()

main()

