#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:     This script is identical to hb.py except that there is no
#              PAIR_TYPES
#
# Author:      user
#
# Created:     27/12/2010
# Copyright:   (c) user 2010
# Licence:     <your licence>
#-------------------------------------------------------------------------------
D_MIN=0
D_MAX=120
D_DELT=2

#R=8.314472
R=0.00198721
T=298


import pickle
import math
import os
import sys
import re

from math import floor
from pmf_atom_typer import LIG_TYPES
from pmf_atom_typer import PRO_TYPES

#from simple_atom_typer import LIG_TYPES
#from simple_atom_typer import PRO_TYPES

def GetDCTypes():
	v=[]
	for i in range(len(LIG_TYPES)):
		for j in range(len(PRO_TYPES)):
			hb_type='%s-%s'%(LIG_TYPES[i],PRO_TYPES[j])
			v.append(hb_type)
	return v

DC_TYPES=GetDCTypes()

def _Smooth(v):
	u=[]
	for d10 in range(D_MIN,D_MAX,D_DELT):
		d=d10*0.1
		k=GetSrIdx(d)
		val=v[k]
		u.append(val)
	w=[]
	for i in range(len(u)):
		val=u[i]
		if i>=8 and i<len(u)-8:
			s=0.0
			for j in range(17):
				s+=u[i-8+j]
			val=s/17
		w.append(val)
	return w

class PTS:
	def __init__(self,dc):
		self.vect=[]   #[[[A,A0,nij,nij_bulk,rho_ij,rho_bulk,fj,i,j],...],...]
		VR=4*math.pi/3.0*((D_MAX*0.1)**3-(D_MIN*0.1)**3)
##		Vd=[]
##		d=4.0
##		for l in range(9):
##			vol=4*math.pi/3.0*d**3
##			if l>0:
##				vol=4*math.pi/3.0*((d+1)**3-d**3)
##				Vd.append(vol)
##				d+=1.0
		#to calculate volume correction factor
		v_fj=[]
		for i in range(len(LIG_TYPES)):
			u=[]
			nkj_bulk=sum(dc.v_nkj[i])
			nlj_bulk=sum(dc.v_nlj[i])
			vj_bulk=0.0
			if (nkj_bulk+nlj_bulk)>0:
				vj_bulk=float(nkj_bulk)/(nkj_bulk+nlj_bulk)
			for l in range(9):
				nkj=dc.v_nkj[i][l]
				nlj=dc.v_nlj[i][l]
				vj=0.0
				if (nkj+nlj) >0:
					vj=float(nkj)/(nkj+nlj)
				fj=1.0
				if vj>1E-5:
					fj=vj_bulk/vj
				u.append(fj)
			w=_Smooth(u)
			v_fj.append(w)
		#nall_bulk=dc.nall_bulk
		#to calculate protein-ligand atom pair potential A_ij(r)
		for j in range(len(DC_TYPES)):
			i=j/len(PRO_TYPES) #lig atom type idx
			nij_bulk=sum(dc.v_nij[j])
			w=[]
			for k in range(len(dc.v_nij[j])):
				d=(D_MIN+k*D_DELT)*0.1
				d1=(D_MIN+(k+1)*D_DELT)*0.1
				A=0
				nij=dc.v_nij[j][k]
				#nall=dc.v_nall[k]
				vol=4*math.pi/3.0*(d1**3-d**3)
				#print d,Vd,VR,nij,nij_bulk
				rho_ij=nij/vol
				rho_bulk=nij_bulk/VR

				fj=v_fj[i][k]
				A0=0
				if nij_bulk<1000:
					A=0
					A0=0
				elif float(nij)/nij_bulk<1E-5:
					A=3
					A0=3
##				elif float(nij)/nij_bulk<3.5E-4:#1E-3 : #Set this range because the distribution of nij is highly biased
##					A=4000
##					A0=4000
				else:
					#if fj<1E-5:fj=1.0
					#print DC_TYPES[j],jj,nij,nij_bulk,Vd,vc.VR,rho_bulk,rho_ij,fj
					A=R*T*(math.log(rho_bulk)-math.log(rho_ij)-math.log(fj))
					A0=R*T*(math.log(rho_bulk)-math.log(rho_ij))
				#print nij,nij_bulk,nall,nall_bulk
				w.append([A,A0,nij,nij_bulk,rho_ij,rho_bulk,fj,i,j])
			self.vect.append(w)
	def __str__(self):
		lines=''
		for j in range(len(DC_TYPES)):
			d_type='%s'%(DC_TYPES[j])
			lines+='A_ij %s\n'%(d_type)
			for k in range(len(self.vect[j])):
				lines+=' %7.1f'%((D_MIN+(k+1)*D_DELT)*0.1)
			lines+='\n'
			for k in range(len(self.vect[j])):
				A=self.vect[j][k][0]
				lines+=' %7.1f'%A
			lines+='\n'
		lines+='\n'
		return lines
##	def Smooth12421(self):
##		v=[]
##		for j in range(len(DC_TYPES)):
##			w=[]
##			for k in range(len(self.vect[j])):
##				A0=self.vect[j][k][0]
##				A_1=A0
##				A_2=A0
##				A1=A0
##				A2=A0
##				if k>=1:
##					A_1=self.vect[j][k-1][0]
##				if k>=2:
##					A_2=self.vect[j][k-2][0]
##				if k<(len(self.vect[j])-1):
##					A1=self.vect[j][k+1][0]
##				if k<(len(self.vect[j])-2):
##					#print k+2
##					A2=self.vect[j][k+2][0]
##				A=(A_2+2*A_1+4*A0+2*A1+A2)*0.1
##				w.append(A)
##			v.append(w)
##		for j in range(len(DC_TYPES)):
##			for k in range(len(self.vect[j])):
##				self.vect[j][k][0]=v[j][k]

def PrintPotentialBatch(pts,dir_name='temp'):
	if not os.path.exists(dir_name):
		os.mkdir(dir_name)
	for j in range(len(DC_TYPES)):
		d_type='%s'%(DC_TYPES[j])
		#print d_type,i,j
		#lines+='A_ij %s\n'%(d_type)
		ap1,ap2=d_type.split('-')
		fo_name='%s/%s_%s.csv'%(dir_name,ap1,ap2)
		if os.path.exists(fo_name):
			print 'File %s exists, check and run the script again.'%(fo_name)
			sys.exit(0)
		fo=open(fo_name,'wb')
		for k in range(len(pts.vect[j])):
			d=(D_MIN+(k+1)*D_DELT)*0.1
			A=pts.vect[j][k][0]
			A0=pts.vect[j][k][1]
			nij=pts.vect[j][k][2]
			nij_bulk=pts.vect[j][k][3]
			rho_ij=pts.vect[j][k][4]
			rho_bulk=pts.vect[j][k][5]
			fj=pts.vect[j][k][6]
			lig_type=LIG_TYPES[pts.vect[j][k][7]]
			dc_type=DC_TYPES[pts.vect[j][k][8]]
			fo.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'\
%(d,A,A0,nij,nij_bulk,rho_ij,rho_bulk,fj,lig_type,dc_type))
		fo.close()
def GetSrIdx(d):
	idx=0
	if d>=4:
		idx=int(d)-3
	return idx

class DC:
	def __init__(self):
		self.v_nij=[]      #length=34*17   numbers of atom pair ij in the sphere shell

		self.v_nkj=[] #length=34*9   numbers of protein atoms k of any type around a ligand atom of type j.  GetStrIdx()
		self.v_nlj=[] #length=34*9   numbers of ligand atoms l of any type around a ligand atom of type j.  GetStrIdx()

		for j in range(len(DC_TYPES)):
			self.v_nij.append([0]*len(range(D_MIN,D_MAX,D_DELT)))

		for j in range(len(LIG_TYPES)):
			self.v_nkj.append([0]*9)
			self.v_nlj.append([0]*9)

	def __str__(self):
		lines=''
		for j in range(len(DC_TYPES)):
			lines+='N_ij %s\n'%(DC_TYPES[j])
			for k in range(len(self.v_nij[j])):
				lines+=' %7.1f'%((D_MIN+(k+1)*D_DELT)*0.1)
			lines+='\n'
			for k in range(len(self.v_nij[j])):
				occ=self.v_nij[j][k]
				lines+=' %7d'%occ
			lines+='\n'
			nij_bulk=sum(self.v_nij[j])
			lines+='%d\n'%nij_bulk
			lines+='N_ij_bulk %s %d\n'%(DC_TYPES[j],nij_bulk)
		lines+='\n'
##		for i in range(len(LIG_TYPES)):
##			lines+='N_kj %s\n'%(LIG_TYPES[i])
##			for k in range(len(self.v_nkj_raw[i])):
##				lines+=' %7.1f'%(k+4)
##			lines+='\n'
##			for k in range(len(self.v_nkj_raw[i])):
##				occ=self.v_nkj_raw[i][k]
##				lines+=' %7d'%occ
##			lines+='\n'
##			lines+='%d\n'%(sum(self.v_nkj_raw[i]))
##			lines+='N_kj_bulk %s %d\n'%(LIG_TYPES[i],self.v_nkj_bulk[i])
##		lines+='\n'
##		for i in range(len(LIG_TYPES)):
##			lines+='N_lj %s\n'%(LIG_TYPES[i])
##			for k in range(len(self.v_nlj_raw[i])):
##				lines+=' %7.1f'%(k+4)
##			lines+='\n'
##			for k in range(len(self.v_nlj_raw[i])):
##				occ=self.v_nlj_raw[i][k]
##				lines+=' %7d'%occ
##			lines+='\n'
##			lines+='%d\n'%(sum(self.v_nlj_raw[i]))
##			lines+='N_lj_bulk %s %d\n'%(LIG_TYPES[i],self.v_nlj_bulk[i])
		lines+='\n\n'
		return lines

	def Add(self, raw_file_name=''):
		raw_lines=_DCFileReader(raw_file_name)
		for (lt_idx,pt_idx,dc_idx,d,IS_LIG_LIG) in raw_lines:
			if d<D_MIN*0.1 or d>=D_MAX*0.1:
				continue
			k=int((d*10.0-D_MIN)/D_DELT)
			l=GetSrIdx(d)
			#print i,j,k,l,d
			if not IS_LIG_LIG:
				self.v_nij[dc_idx][k]+=1
				if pt_idx!=16: #only count non-hydrogen protein atoms
					self.v_nkj[lt_idx][l]+=1
			else:
				self.v_nlj[lt_idx][l]+=1
		
	def Add2(self, dc):
		for j in range(len(DC_TYPES)):
			for k in range(len(self.v_nij[j])):
				self.v_nij[j][k]+=dc.v_nij[j][k]
		for i in range(len(LIG_TYPES)):
			for k in range(len(self.v_nkj[i])):
				self.v_nkj[i][k]+=dc.v_nkj[i][k]
				self.v_nlj[i][k]+=dc.v_nlj[i][k]

def UnformatedInput(fn):
	f=open(fn)
	dc=pickle.load(f)
	f.close()
	return dc

def UnformatedOutput(dc,fn):
	f=open(fn,'wb')
	pickle.dump(dc,f)
	f.close()

def _DCParseLine(line):
	IS_LIG_LIG=False
	v=line.rstrip().split(',')
	dt_type=v[0]
#	ids=re.findall('<(\d+)>',v[1])
	dt_idx=-1

	lt_idx=eval(v[1])
	pt_idx=eval(v[2])
	if lt_idx<0:
		print 'Ligand type is not found. The following line has been neglected:\n%s'%(line)
		return None
	if pt_idx<0:
		print 'DC type %s is not found. The followint line has been neglected:\n%s'%(dt_type,line)
		return None
	if dt_type[0]=='#':
		IS_LIG_LIG=True
		dt_type=dt_type[1:]
	else:
		dt_idx=DC_TYPES.index(dt_type)
	d=eval(v[4][:-1])
	return (lt_idx,pt_idx,dt_idx,d,IS_LIG_LIG)

def _DCGetLines(w):
	v=[]
	for line in w:
		#if line[0]=='#':continue
		data=_DCParseLine(line)
		if data:
			v.append(_DCParseLine(line))
	return v

def _DCFileReader(fn=''):
	w=open(fn).readlines()
	v=_DCGetLines(w)
	return v


def main():

	dc0=DC()
	#print dc0


	dc0.Add('test_raw_print_dist.out')
	#print dc0

	pts1=PTS(dc0)
	print pts1

if __name__ == '__main__':
    main()
