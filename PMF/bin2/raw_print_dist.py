#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      user
#
# Created:     23/12/2010
# Copyright:   (c) user 2010
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import pybel
import openbabel as ob
import math
import pmf_atom_typer as pat
#import simple_atom_typer as pat

ob.obErrorLog.SetOutputLevel(0)

def GetDistance(a1,a2):
##a1=donor H
##a2=acceptor
	dis=a1.GetDistance(a2)
	return dis


def FormatedOutput(atom1,atom2):
	res=atom2.GetResidue()
	res_name='UNK'
	res_num=-1
	if res:
		res_name=res.GetName()
		res_num=res.GetNum()
	line="[LIG:%s<%d>...%s%d:%s<%d>]"%(atom1.GetType(),atom1.GetIdx()\
,res_name,res_num,atom2.GetType(),atom2.GetIdx())
	return line

def RawDisPrint(lig,pro):
	v=[]
	for a1 in lig.atoms:
		atom1=a1.OBAtom
##		if a1.OBAtom.IsHydrogen():
##			continue
		a1_type=pat.LigAtomTyper(atom1)
		a1_idx=-1
		if a1_type in pat.LIG_TYPES:
			a1_idx=pat.LIG_TYPES.index(a1_type)
		for a2 in pro.atoms:
			atom2=a2.OBAtom
##			if a2.OBAtom.IsHydrogen():
##				continue
			distance=GetDistance(atom1,atom2)
			if distance>=12.0:continue
			a2_type=pat.ProAtomTyper(atom2)
			a2_idx=-1
			if a2_type in pat.PRO_TYPES:
				a2_idx=pat.PRO_TYPES.index(a2_type)
			dis_type='%s-%s'%(a1_type,a2_type)
			info=FormatedOutput(atom1,atom2)
			line='%s,%d,%d,%s,%.4fA\n'%(dis_type,a1_idx,a2_idx,info,distance)
			v.append(line)
		for a2 in lig.atoms:
			atom2=a2.OBAtom
##			if a2.OBAtom.IsHydrogen():
##				continue
			if atom1.GetIdx()==atom2.GetIdx():
				continue
			distance=GetDistance(atom1,atom2)
			if distance>=12.0:continue
			a2_type=pat.LigAtomTyper(atom2)
			a2_idx=-1
			if a2_type in pat.LIG_TYPES:
				a2_idx=pat.LIG_TYPES.index(a2_type)
			dis_type='#%s-%s'%(a1_type,a2_type)
			info=FormatedOutput(atom1,atom2)
			line='%s,%d,%d,%s,%.4fA\n'%(dis_type,a1_idx,a2_idx,info,distance)
			v.append(line)
	return v

def RawPrint(lig_fn,pro_fn,output_fn):
	lig=pybel.readfile('sdf',lig_fn).next()
	pro=pybel.readfile('pdb',pro_fn).next()
	fo=open(output_fn,'wb')
	v=RawDisPrint(lig,pro)
	fo.writelines(v)
	fo.close()

def main():
	lig_fn='1AQI_lig.sdf'
	pro_fn='1AQI_pro.pdb'
	output_fn='1AQI_raw_print_dist.out'
	RawPrint(lig_fn,pro_fn,output_fn)


if __name__ == '__main__':
	main()
