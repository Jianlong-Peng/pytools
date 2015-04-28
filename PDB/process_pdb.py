'''
#=============================================================================
#     FileName: process_pdb.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-06-07 22:02:03
#   LastChange: 2014-06-25 00:43:20
#      History:
#=============================================================================
'''
import sys
import tools

def main(argv=sys.argv):
    if len(argv) != 4:
        print "\n  Usage: %s [option] out_dir"%argv[0]
        print "\n  [option]"
        print "    -f pdb_file"
        print "    -l file.  pdb file in each line"
        print "  out_dir: where to save extracted protein(s), ligand(s), and water(s)"
        print "\n  a file named 'names_pro_lig.log' will be generated"
        print ""
        sys.exit(1)

    
    outf_log = open("names_pro_lig.log","a")
    if argv[1] == "-f":
        tools.process(argv[2], argv[3], outf_log)
    elif argv[1] == "-l":
        inf = open(argv[2],"r")
        for line in inf:
            if line.startswith('#'):
                continue
            tools.process(line.strip(), argv[3], outf_log)
        inf.close()
    else:
        print "Error: invalid option ",argv[1]
    outf_log.close()

main()

