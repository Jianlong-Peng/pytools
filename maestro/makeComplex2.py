'''
#=============================================================================
#     FileName: makeComplex.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-06-21 15:12:17
#   LastChange: 2013-06-24 09:32:03
#      History:
#=============================================================================
'''

import sys
import os
try:
    from schrodinger import structure, structureutil
    import stageTools
except ImportError:
    print "\n Please run `%s\\maestro_variable.bat` first!!!\n"%(sys.argv[0][:sys.argv[0].rfind("\\")])
    sys.exit(1)

receptor = {}
ligands = []

def parse_config(config_file):
    inf = open(config_file,"r")
    global receptor
    global ligands
    line = inf.readline()
    while line != "":
        if line.startswith('#') or line.strip()=="":
            line = inf.readline()
            continue
        if line.startswith("receptor"):
            if receptor:
                print "  Error: more than one receptor is given!!!!"
                inf.close()
                sys.exit(1)
            while line!="" and line.strip()!="":
                line = line.split()
                receptor[line[0][:-1]] = line[1]
                line = inf.readline()
        elif line.startswith("ligand"):
            temp = {}
            while line!="" and line.strip()!="":
                line = line.split()
                temp[line[0][:-1]] = line[1]
                line = inf.readline()
            ligands.append(temp)
        else:
            print r"  Error: invalid line `%s`"%line
            inf.close()
            sys.exit(1)
    inf.close()

def writeInpFile(outfile, complex_file):
    '''
    outfile: xxx.inp, used to run covalent docking
    complex_file: while complex(generated here) locates
    '''
    outf = open(outfile,"w")
    print >>outf, """

[ SET:COMPLEXES]
    VARCLASS Structures
    FILES %s

[ STAGE:PRIME ]
    STAGECLASS              prime.PrimeStage
    INPUTS                  COMPLEXES
    OUTPUTS                 DOCKED_OUT
    PRIME_TYPE              COVALENT_DOCKING
    LIGAND                  X:1
    INCLUDE_RESIDUE         yes
    NUM_OUTPUT_STRUCT       1


[ USEROUTS ]
    USEROUTS DOCKED_OUT
    STRUCTOUT DOCKED_OUT"""%complex_file
    
    outf.close()


def main(argv=sys.argv):
    if len(argv) != 2:
        print "\n    OBJ: to make complex between given ligand and receptor"
        print "\n  Usage: makeComplex.py config.txt"
        print "  config.txt: where information about receptor and ligand locate"
        sys.exit(1)

    global receptor
    global ligands
    print "to parse `%s`"%argv[1]
    parse_config(argv[1])

    try:
        receptor_file = receptor['receptor_file']
        receptor_leaving_atom = int(receptor['receptor_leaving_atom'])
        receptor_staying_atom = int(receptor['receptor_staying_atom'])
    except KeyError:
        print "Error: missing keys for receptor!"
        sys.exit(1)
    receptor_st = structure.StructureReader(receptor_file).next()
    
    print "receptor: `%s`"%receptor_file
    print "receptor_leaving_atom: `%d`"%receptor_leaving_atom
    print "receptor_staying_atom: `%d`\n"%receptor_staying_atom

    for i,ligand in enumerate(ligands):
        print "%d. to process ligand(s) from file %s"%(i+1,ligand['ligand_file'])
        try:
            ligand_file = ligand['ligand_file']
            smarts = ligand['smarts']
            complex_out = ligand['complex_out']
            #base_name,ext = os.path.splitext(complex_out)
        except KeyError:
            print "  Error: missing key(s)!"
        else:
            if os.path.exists(complex_out):
                print "  Error: `%s` is already exists!"%complex_out
                print "         The corresponding complexes will not be generated!!"
                continue
            total_success = 0
            outf = structure.StructureWriter(complex_out)
            inf = structure.StructureReader(ligand_file)
            for j,ligand_st in enumerate(inf):
                if ligand_st.title == '':
                    print "  > process %dth ligand"%(j+1)
                else:
                    print "  > process %dth ligand (%s)"%(j+1,ligand_st.title)
                match_list = structureutil.evaluate_smarts(ligand_st,smarts)
                print "    totally found %d matches"%len(match_list)
                for matched_atoms in match_list:
                    ligand_reactive_atom = matched_atoms[0]
                    print "    - try to make bond between atom `%d`(ligand) and `%d`(receptor)"%(ligand_reactive_atom,receptor_staying_atom)
                    complexes_st = stageTools.makeComplex(receptor_st,receptor_leaving_atom,receptor_staying_atom,ligand_st,ligand_reactive_atom)
                    if complexes_st is None:
                        print "      Fail"
                    else:
                        #modification. 2013-06-24 20:19
                        tmp_count = 0
                        for complex_st in complexes_st:
                            outf.append(complex_st)
                            total_success += 1
                            tmp_count += 1
                        print "      Success (%d)"%tmp_count
                        #out_name = "%s_%s_%d%s"%(base_name,title,ligand_reactive_atom,ext)
                        #complex_st.write(out_name)
                        #inp_name = "%s_%s_%d.inp"%(base_name,title,ligand_active_atom)
                        #writeInpFile(inp_name,out_name)
                        #print "      Success, complex has been written to file"
                        #print "        `%s`"%out_name
                        #print "      And, corresponding inp file was generated"
                        #print "        `%s`"%inp_name
            print "  => totally %d complexes have been written to file"%total_success
            print "     `%s`"%complex_out
            base_name,ext = os.path.splitext(complex_out)
            inp_name = base_name+".inp"
            writeInpFile(inp_name,complex_out)
            print "     corresponding inp file is saved in"
            print "     `%s`"%inp_name

main()


