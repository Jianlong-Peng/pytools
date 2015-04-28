#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      BreaShare
#
# Created:     21/05/2013
# Copyright:   (c) user 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

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
    return list

def GetFinger(lig,pro,pts,pdbid):
    mscore=0.0  #with volumn correction
    mscorev=0.0 #without volumn correction
    cutoff=12.0
    Finger=[]
    Fingerv=[]
    f1=open("%s_GlideCfinger.txt"%pdbid,'wb')
    f2=open("%s_GlideCfingerv.txt"%pdbid,'wb')
    for i in range(len(dist.DC_TYPES)):
        Finger.append([0]*4)
    for i in range(len(dist.DC_TYPES)):
        Fingerv.append([0]*4)
    for a1 in lig.atoms:
        if a1.OBAtom.IsHydrogen():
            continue
        a1_type=pat.LigAtomTyper(a1.OBAtom)
        for a2 in pro.atoms:
            if a2.OBAtom.IsHydrogen():
                continue
            a2_type=pat.ProAtomTyper(a2.OBAtom)
            if a2_type=='OW':continue
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
            ascore=pts.vect[j][k][0]
            ascorev=pts.vect[j][k][1]
            m=int(d/3.0)
            Finger[j][m]+=ascore
            Fingerv[j][m]+=ascorev
            mscore+=ascore
            mscorev+=ascorev
    CFinger=CList(Finger)
    CFingerv=CList(Fingerv)
    print CFinger
    for j in range(len(dist.DC_TYPES)):
        dist_type=dist.DC_TYPES[j]
        for i in range(4):
            afinger=CFinger[j][i]
            afingerv=CFingerv[j][i]
            f1.write("%s-%d%7.2f\n"%(dist_type,i,afinger))
            f2.write("%s-%d%7.2f\n"%(dist_type,i,afingerv))
    f1.close()
    f2.close()
    return mscore,mscorev

def GetFileList():
    f=open("toppose.list")
    v=[]
    for line in f.readlines():
        ligname=line.rstrip().lstrip()
       # proname=ligname.split("_")[0]
        lig_fn='%s.sdf'%ligname
        pro_fn='3VO3_glide.pdb'
        pdbid=ligname
        v.append((lig_fn,pro_fn,pdbid))
    return v

def main():
    file_list=GetFileList()
    dc_fn='dist_splited_train_0525.dat'
    dc=dist.UnformatedInput(dc_fn)
    pts=dist.PTS(dc)
    for i in range(len(file_list)):
        lig_fn=file_list[i][0]
        pro_fn=file_list[i][1]
        pdbid=file_list[i][2]
        if os.path.isfile(lig_fn):
            lig=pybel.readfile('sdf',lig_fn).next()
            pro=pybel.readfile('pdb',pro_fn).next()
            score,scorev=GetFinger(lig,pro,pts,pdbid)
            print score,scorev
if __name__=='__main__':
    main()
    
    
    



























































    
