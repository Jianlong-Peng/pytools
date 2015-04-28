import sys
import matplotlib.pyplot as plt
import numpy as np


def main(argv=sys.argv):
    if len(argv) != 4:
        print "\n  Usage: %s active_titles.txt ranking.txt roc\n"%argv[0]
        print "  active_titles.txt: one title per each line"
        print "                     those will be cheated as active ligands!"
        print "  ranking.txt: from 2nd line one, each should be: title docking_score"
        print "               ligands will be sorted according to docking_score!!!"
        print "  roc        : figure"
        sys.exit(1)
    
    inf_active = open(argv[1],"r")
    actives = {}  #title:ranking
    for line in inf_active:
        actives[line.strip()] = 0
    inf_active.close()
    
    inf = open(argv[2],"r")
    line = inf.readline()
    title_score = []
    for line in inf:
        line = line.split()
        title_score.append((line[0],float(line[1])))
    inf.close()
    #title_score.sort(key=lambda x:x[1],reverse=True)  #in case that score has positive values
    title_score.sort(key=lambda x:x[1])   #in case that score has negative values
    num_active = 0
    num_decoy = 0
    i = 0
    tp = 0
    fp = 0
    tps = []
    fps = []  #1-specificity = fp / (tn+fp)
    is_active = []
    for title,score in title_score:
        i += 1
        if actives.has_key(title):
            num_active += 1
            actives[title] = i
            tp += 1
            tps.append(tp)
            fps.append(fp)
            is_active.append(True)
        else:
            num_decoy += 1
            fp += 1
            is_active.append(False)
    
    x = [0.]
    y = [0.]
    x.extend([1.0*item/num_decoy for item in fps])
    y.extend([1.0*item/num_active for item in tps])
    x.append(1.)
    y.append(1.)
    
    #ranking
    title_rank = actives.items()
    title_rank.sort(key=lambda x:x[1])
    print "\n  rank ligand_title"
    for title,rank in title_rank:
        print "  %-4d %-s"%(rank,title)
    
    #AUC
    sum_rel_ranks = 0.
    for title,rank in actives.iteritems():
        sum_rel_ranks += (1.0 * rank / (num_active+num_decoy))
    auac = 1.0 - ((1.0/num_active) * sum_rel_ranks)
    Ri = 1.0*num_decoy / (num_active + num_decoy)
    Ra = 1.0*num_active / (num_active + num_decoy)
    roc = (auac/Ri) - (Ra/(2*Ri))
    print "\n  auac=%g, roc=%g"%(auac,roc)
    
    #enrichment factor
    total_ligands = num_active + num_decoy
    temp = [1,2,5,10,20,30,40,50,60,70,80,90,100]
    print "\n  %samples  EF"
    for value in temp:
        num_select_total = int(1.0 * total_ligands * value / 100)
        num_select_active = is_active[:num_select_total].count(True)
        ef = (1.0 * num_select_active / num_active) / (1.0 * num_select_total / total_ligands)
        print "  %-2d        %g"%(value,ef)
    
    #plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x,y,marker="o",markersize=3)
    
    x = np.linspace(0,1,10)
    y = np.linspace(0,1,10)
    ax.plot(x,y,'-',color='black')
    
    ax.grid()
    ax.set_xlabel("1-specificity",fontname="Times New Roman")
    ax.set_ylabel("sensitivity",fontname="Times New Roman")
    ax.set_xlim(-0.02,1.02)
    ax.set_ylim(-0.02,1.02)
    ax.set_title("AUC=%g"%roc,fontname="Times New Roman")
    
    #plt.show()
    _format = argv[3].split(".")[-1]
    #plt.savefig(argv[3],dpi=120,format=_format)
    plt.savefig(argv[3],format=_format)

    
main()
    





