'''
#=============================================================================
#     FileName: ifSame.py
#         Desc: to judge if two files have same contents
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-08-10 17:29:38
#   LastChange: 2012-08-10 17:44:11
#      History:
#=============================================================================
'''

import sys

def main(argv=sys.argv):
    if len(argv) != 3:
        print "Usage: ifSame.py file1 file2"
        print "blank lines at the begin and end of the file will be ignored!"
        sys.exit(1)

    inf1 = open(argv[1],"r")
    inf2 = open(argv[2],"r")

    line1 = inf1.readline()
    line2 = inf2.readline()

    while line1!="" and line1.strip()=="":
        line1 = inf1.readline()
    while line2!="" and line2.strip()=="":
        line2 = inf2.readline()

    line_no = 1

    while line1!="" and line2!="":
        if line1.strip()!=line2.strip():
            break
        line1 = inf1.readline()
        line2 = inf2.readline()
        line_no += 1

    while line1!="" and line1.strip()=="":
        line1 = inf1.readline()
    while line2!="" and line2.strip()=="":
        line2 = inf2.readline()

    if line1!="" or line2!="":
        print >>sys.stderr, "the two files contain different contents!"
        print >>sys.stderr, "be different begin with line",line_no
    else:
        print "the two files contain same contents!"

    inf1.close()
    inf2.close()

main()



