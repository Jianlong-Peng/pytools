'''
#=============================================================================
#     FileName: plot_num_atom_pairs.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-07-21 14:22:01
#   LastChange: 2014-07-21 16:34:27
#      History:
#=============================================================================
'''
import sys
import dist
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle
from math import log10

lig_types = ["CF","CP","AF","AP","C3","CW","CO","CN","C0","NC","NP","NA","ND","NR","N0","NS","OC","OA","OE","OS","OD","P","SA","SD","SO","HL","F","CL","Br","I"]
pro_types = ["CF","CP","AF","AP","CO","CN","NA","NC","ND","OC","OA","OD","OW","SA","SD","HH","ME"]
def main(argv=sys.argv):
    if len(argv) != 3:
        print "\n  Usage: %s dist_splited_train_0525.dat out.tif"%argv[0]
        print "  out.tif: where to save the image"
        print ""
        sys.exit(1)

    inf = open(argv[1])
    ap = pickle.load(inf)
    inf.close()

    m = [[0 for j in xrange(len(lig_types))] for i in xrange(len(pro_types))]
    for i in xrange(len(pro_types)):
        for j in xrange(len(lig_types)):
            types = '%s-%s'%(lig_types[j], pro_types[i])
            try:
                index = dist.DC_TYPES.index(types)
            except ValueError:
                print types
                sys.exit(1)
            try:
                m[i][j] = log10(sum(ap.v_nij[index]))
            except ValueError:
                m[i][j] = 0.

    m = np.asarray(m)

    matplotlib.rcParams['font.serif'] = 'Times New Roman'
    matplotlib.rcParams['figure.figsize'] = (12,6)
    #plt.pcolormesh(m)
    heatmap = plt.pcolor(m)
    plt.xlim(0, len(lig_types))
    plt.ylim(0, len(pro_types))
    plt.xticks(np.arange(m.shape[1])+0.5, lig_types, fontsize='x-small')
    plt.xlabel("Ligand Atom Types", fontsize='small')
    plt.yticks(np.arange(m.shape[0])+0.5, pro_types, fontsize='x-small')
    plt.ylabel("Protein Atom Types", fontsize='small')
    plt.grid()
    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    cb = plt.colorbar()
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.pcolor(m)
    ax.grid()
    #ax.set_frame_on(True)
    ax.set_xlim(0, len(atom_type.LTYPE[0]))
    ax.set_xticks(np.arange(m.shape[1])+0.5, minor=False)
    ax.set_xticklabels(atom_type.LTYPE[0], minor=False, fontsize='xx-small')
    ax.set_xlabel('Ligand Atom Types', fontsize='small')
    ax.set_ylim(0, len(atom_type.PTYPE[0]))
    ax.set_yticks(np.arange(m.shape[0])+0.5, minor=False)
    ax.set_yticklabels(atom_type.PTYPE[0], fontsize='x-small')
    ax.set_ylabel('Protein Atom Types', fontsize='small')

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    """


    plt.savefig(argv[2],format='tif',dpi=150)

main()

