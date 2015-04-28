'''
#=============================================================================
#     FileName: splitSdf.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-06-01 21:01:44
#   LastChange: 2013-06-01 21:18:45
#      History:
#=============================================================================
'''
import sys
import os

def main(argv=sys.argv):
    if len(argv)!=3 and len(argv)!=4:
        print "\n   OBJ: to split 'in.sdf' into pieces and saved {in out_dir}"
        print "        each file with one molecule!"
        print "\n Usage: splitSdf.py in.sdf out_dir [name]\n"
        print " [name]: specify the content used as output file name"
        print "   if not given, content in the 1st line will be used as file name!"
        print "   otherwise, the property{name} in '> <...>' will be used as file name!"
        print " if can't find {name} or it's empty, then kind of 'molecule_X.sdf' will be used\n"
        sys.exit(1)

    assert os.path.isfile(argv[1])
    if not os.path.exists(argv[2]):
        os.mkdir(argv[2])
    inf = open(argv[1],"r")
    line = inf.readline()
    count = 0
    while line != "":
        count += 1
        name = ""
        if len(argv) == 3:
            name = line.strip()
        info = line
        line = inf.readline()
        while line.strip() != "$$$$":
            info += line
            if len(argv) == 4:
                if len(line)>3 and line[:3]=="> <" and line[3:-2]==argv[3]:
                    line = inf.readline()
                    info += line
                    name = line
            line = inf.readline()
        info += line
        if name == "":
            name = "molecule_%d"%count
        outf = open(os.path.join(argv[2],name+".sdf"),"w")
        outf.write(info)
        outf.close()
        #go to next molecule
        line = inf.readline()
    inf.close()
    print "Totally %d molecules read and saved!"%count

main()

