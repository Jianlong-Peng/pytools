'''
#=============================================================================
#     FileName: rmDesFromLibsvm.py
#         Desc: to remove descriptor(s) from a libsvm file
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-09-14 13:43:22
#   LastChange: 2012-09-14 13:56:19
#      History:
#=============================================================================
'''

import sys

def main(argv=sys.argv):
    if len(argv) < 4:
        print "Usage: rmDesFromLibsvm.py in.svm out.svm tag1[,tag2,...]"
        sys.exit(1)

    rm_tags = set([int(_) for _ in argv[3:]])
    inf = open(argv[1],"r")
    outf = open(argv[2],"w")
    num_tags = 0
    for line in inf:
        tmp = int(line.split()[-1].split(":")[0])
        if tmp > num_tags:
            num_tags = tmp
    tags_map = [0 for i in xrange(num_tags)]
    inf.seek(0)
    for line in inf:
        line = line.split()
        outf.write(line[0])
        i = 1
        j = 0
        for item in line[1:]:
            j += 1
            item = item.split(":")
            if int(item[0]) not in rm_tags:
                outf.write(" "+str(i)+":"+item[1])
                tags_map[j-1] = i
                i += 1
        outf.write("\n")

    inf.close()
    outf.close()

    print "orig_tag -> now_tag"
    for i in xrange(len(tags_map)):
        if tags_map[i]:
            print "%3d -> %-d"%(i+1,tags_map[i])
        else:
            print "%3d -> --" %(i+1)

main()



