'''
#=============================================================================
#     FileName: addDes2Libsvm.py
#         Desc: add one more descriptors to an existing libsvm file
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-08-13 11:13:45
#   LastChange: 2012-10-22 14:21:45
#      History:
#=============================================================================
'''

import sys

def main(argv=sys.argv):
    if len(argv)!=4 and len(argv)!=5:
        print "Usage addDes2Libsvm.py libsvm_file des_file out_libsvm [begin_with_col]"
        print "libsvm_file: an existing libsvm file"
        print "des_file:"
        print "  1st line: titles"
        print "  2~* line: des1 des2 ..."
        print "Suppose that mols in `des_file` should be correspond to mols `libsvm_file`"
        print "[begin_with_col]: <default: 2>"
        print "  the first {begin_with_col-1} values are ignored in each line."
        sys.exit(1)

    if len(argv) == 5:
        begin_col = int(argv[4])
    else:
        begin_col = 2

    inf_des = open(argv[2],"r")
    des = []
    line = inf_des.readline()
    for line in inf_des:
        des.append(line.split()[begin_col-1:])
    inf_des.close()

    inf_libsvm = open(argv[1],"r")
    num_features = 0
    num_instances = 0
    for line in inf_libsvm:
        line = line.split()
        max_features = int(line[-1].split(":")[0])
        if num_features < max_features:
            num_features = max_features
        num_instances += 1

    if len(des) != num_instances:
        print >>sys.err, "There are different number of instances in",argv[1],"and",argv[2]
        inf_libsvm.close()
        sys.exit(1)
    
    inf_libsvm.seek(0)
    outf = open(argv[3],"w")
    i = 0
    for line in inf_libsvm:
        outf.write(line.strip())
        index = num_features+1
        for new_des in des[i]:
            outf.write(" %d:%s"%(index,new_des))
            index += 1
        outf.write("\n");
        i += 1

    inf_libsvm.close()
    outf.close()

main()


