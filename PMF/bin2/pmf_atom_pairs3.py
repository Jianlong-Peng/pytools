#!/usr/bin/env python
'''
#=============================================================================
#     FileName: pmf_atom_pairs3.py
#         Desc: based on 'gen_pmf_fingerA.py'
#       Author: jlpeng
#        Email: 
#     HomePage: 
#      Version: 0.0.1
#   LastChange: 2014-08-19 13:06:03
#      History:
#=============================================================================
'''
import sys
import os
import pybel
import dist
import pmf_atom_typer as pat
import raw_print_dist as rpd
from math import floor

def CList(list):
    for each in list:
        each[1]=each[0]+each[1]
        each[2]=each[1]+each[2]
        each[3]=each[2]+each[3]
        each[4]=each[3]+each[4]
        each[5]=each[4]+each[5]
        each[6]=each[5]+each[6]
        each[7]=each[6]+each[7]
        each[8]=each[7]+each[8]
        each[9]=each[8]+each[9]
        each[10]=each[9]+each[10]
        each[11]=each[10]+each[11]
    return list

def GetFinger(lig,pro,pts,pdbid):
    mscore=0.0  #without volumn correction
    mscorev=0.0 #with volumn correction
    cutoff=12.0
    Finger=[]  #without volumn correction
    Fingerv=[]  #with volumn correction
    #f1=open("%s_GlideAfinger.txt"%pdbid,'wb')
    f2=open("%s.num"%pdbid,'wb')
    #for i in range(len(dist.DC_TYPES)):
    #    Finger.append([0]*12)
    #for i in range(len(dist.DC_TYPES)):
    #    Fingerv.append([0]*12)
    for a1 in lig.atoms:
        if a1.OBAtom.IsHydrogen():
            continue
        a1_type=pat.LigAtomTyper(a1.OBAtom)
        if a1_type == '': continue
        for a2 in pro.atoms:
            if a2.OBAtom.IsHydrogen():
                continue
            a2_type=pat.ProAtomTyper(a2.OBAtom)
            if a2_type in ('OW',''):continue
            dis_type='%s-%s'%(a1_type,a2_type)
            if dis_type not in dist.DC_TYPES:
                print 'Atom pair type %s has not been defined, the following pair has been neglected.'%(dis_type)
                print a1.OBAtom.GetType(),a1.OBAtom.GetIdx(),a2.OBAtom.GetType(),a2.OBAtom.GetIdx()
                continue
            j=dist.DC_TYPES.index(dis_type)
            d=rpd.GetDistance(a1.OBAtom,a2.OBAtom)
            if d>=cutoff:
                continue
            k=int((d-dist.D_MIN*0.1)/dist.D_DELT*10.0)
            #print >>f2, "%s,%d,%.6f,%d"%(dis_type,j,d,k)
            print >>f2, "%s,%d,%s,%d,%.6f"%(a1_type,a1.OBAtom.GetIdx(),a2_type,a2.OBAtom.GetIdx(),d)
            mscorev += pts.vect[j][k][0]
            mscore += pts.vect[j][k][1]
    '''
            ascore=pts.vect[j][k][0]
            ascorev=pts.vect[j][k][1]
            m=int(d/1.0)
            Finger[j][m]+=ascore
            Fingerv[j][m]+=ascorev
            mscore+=ascore
            mscorev+=ascorev
    CFinger=CList(Finger)
    CFingerv=CList(Fingerv)
    print CFinger
    for j in range(len(dist.DC_TYPES)):
        dist_type=dist.DC_TYPES[j]
        i=1
        while i<12:
            afinger=CFinger[j][i]
            afingerv=CFingerv[j][i]
            f1.write("%s-%d%10.2f\n"%(dist_type,i-1,afinger))
            f2.write("%s-%d%10.2f\n"%(dist_type,i-1,afingerv))
            i=i+1
    f1.close()
    '''
    f2.close()
    return mscore,mscorev

def GetFileList(infile):
    #f=open("toppose.list")
    f = open(infile,'r')
    v=[]
    for line in f.readlines():
        proname,ligname = line.split()
        pdbid = os.path.basename(ligname).split("_")[0]
        v.append((ligname, proname, pdbid))
    f.close()
    return v

def main(argv=sys.argv):
    if len(argv) != 2:
        print "\n  Usage: %s in.list"%argv[0]
        print "  in.list: each line should be 'proname ligname'"
        print ""
        sys.exit(1)
    file_list=GetFileList(argv[1])
    #to get atom pair potentional A_ij(r)
    dirname = os.path.dirname(argv[0])
    dc_fn = os.path.join(dirname,'dist_splited_train_0525.dat')
    dc=dist.UnformatedInput(dc_fn)
    pts=dist.PTS(dc)
    for i in range(len(file_list)):
        lig_fn=file_list[i][0]
        pro_fn=file_list[i][1]
        pdbid=file_list[i][2]
        lig_format = lig_fn[lig_fn.rfind(".")+1:]
        pro_format = pro_fn[pro_fn.rfind(".")+1:]
        if os.path.isfile(lig_fn):
            lig=pybel.readfile(lig_format,lig_fn).next()
            pro=pybel.readfile(pro_format,pro_fn).next()
            score,scorev = GetFinger(lig,pro,pts,pdbid)
            print pdbid,score,scorev

if __name__=='__main__':
    main()
    

