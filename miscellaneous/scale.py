'''
#=============================================================================
#     FileName: scale.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-12-05 09:44:33
#   LastChange: 2014-02-14 14:20:34
#      History:
#=============================================================================
'''
import sys
from getopt import getopt

def read_parameters(infile):
    min_max = [[],[]]
    inf = open(infile,"r")
    line = inf.readline()
    for line in inf:
        _min,_max = map(float,line.split())
        min_max[0].append(_min)
        min_max[1].append(_max)
    inf.close()
    return min_max

def calc_parameters(infile,sep):
    inf = open(infile,"r")
    line = inf.readline()
    if sep == "":
        num_desc = len(line.split())-1
    elif sep == "tab":
        num_desc = len(line.split("\t"))-1
    else:
        num_desc = len(line.strip().split(sep))-1
    min_max = [[1e38 for _ in xrange(num_desc)],[-1e38 for _ in xrange(num_desc)]]
    for line in inf:
        if sep == "":
            line = line.split()
        elif sep == "tab":
            line = line.strip().split("\t")
        else:
            line = line.strip().split(sep)
        for i in xrange(1,len(line)):
            val = float(line[i])
            if val < min_max[0][i-1]:
                min_max[0][i-1] = val
            if val > min_max[1][i-1]:
                min_max[1][i-1] = val
    inf.close()
    return min_max

def save_parameters(min_max,outfile):
    outf = open(outfile,"w")
    print >>outf, "min\tmax"
    for i in xrange(len(min_max[0])):
        print >>outf, "%.6g\t%.6g"%(min_max[0][i],min_max[1][i])
    outf.close()

def scale(desc, _min, _max):
    if _min == _max:
        return 1.
    return (desc - _min) / (_max - _min)

def main(argv=sys.argv):
    if len(argv) not in (5,7):
        print "\n  Usage: %s [option] input output\n"%argv[0]
        print "  [option]"
        print "  --sep delimiter: the character used to seperate each item"
        print "    tab to specify the tab character"
        print "    <default: space>"
        print "  -s save_filename: save scaling parameters"
        print "  -r restore_filename: restore scaling parameters\n\n"
        print "  . scaling method: (X-X_min)/(X_max-X_min)"
        print "  . if `-s` is given, `infile` will be scalled to (-1,1),"
        print "    and parameters will be saved in `save_filename`"
        print "  . if `-r` is given, scaling `infile` using `restore_filename` instead.\n"
        sys.exit(1)
    
    sep = ""
    save_para = ""
    restore_para = ""

    opts,args = getopt(argv[1:],'s:r:',["sep="])
    if len(args) != 2:
        print "Error: invalid number of arguments ",args
        sys.exit(1)
    for opt,val in opts:
        if opt == "--sep":
            sep = val
        elif opt == "-s":
            save_para = val
        elif opt == "-r":
            restore_para = val
        else:
            print "Error: invalid option",opt
            sys.exit(1)

    if save_para=="" and restore_para=="":
        print "Error: one of -s,-r is needed!!"
        sys.exit(1)
    if save_para!="" and restore_para!="":
        print "Error: only one of -s,-r is allowed!!"
        sys.exit(1)

    min_max = None  #[[min1,min2,...],[max1,max2,...]]
    if restore_para != "":
        min_max = read_parameters(restore_para)
    if save_para != "":
        min_max = calc_parameters(args[0],sep)
        save_parameters(min_max,save_para)

    inf = open(args[0],"r")
    line = inf.readline()
    if sep == "":
        num_desc = len(line.split())-1
    elif sep == "tab":
        num_desc = len(line.split("\t"))-1
    else:
        num_desc = len(line.strip().split(sep))-1
    if(num_desc != len(min_max[0])):
        print "Error: incompatible number of descriptors from",restore_para
        inf.close()
        sys.exit(1)
    outf = open(args[1],"w")
    outf.write(line)
    for line in inf:
        if sep == "":
            line = line.split()
        elif sep == "tab":
            line = line.strip().split("\t")
        else:
            line = line.strip().split(sep)
        outf.write(line[0])
        desc = map(float,line[1:])
        for i in xrange(len(desc)):
            val = scale(desc[i],min_max[0][i],min_max[1][i])
            if sep == "":
                outf.write("\t%.6g"%val)
            elif sep == "tab":
                outf.write("\t%.6g"%val)
            else:
                outf.write("%s%.6g"%(sep,val))
        outf.write("\n")
    inf.close()
    outf.close()

main()


