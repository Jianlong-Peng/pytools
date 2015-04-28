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
ob.obErrorLog.SetOutputLevel(0)

LIG_TYPES=['CF','CP','AF','AP','C3','CW','CO','CN','C0','NC','NP','NA','ND','NR','N0','NS','OC','OA','OE','SO','OS', 'OD','P','SA','SD','HL','CL','F', 'Br','I']
PRO_TYPES=['CF','CP','AF','AP','CO','CN','NC','ND','ME','OC','OA','OD','OW','SA','SD','NA','HH']

##-----
##Functions: get number of hydrogen atoms connected to an atom
##Input: the Ob Object: OBAtom()
##Output: none
##Return: an int, number of hydrogens
##Calls: none
##Called By: ProTyper(); LigTyper()
def NumConnectH(a):
	num_connected_h=0
	for an in ob.OBAtomAtomIter(a):
		if an.IsHydrogen():
			num_connected_h+=1
	return num_connected_h

##----
##Functions: decide a nitrogen atom is a charged nitrogen or not
##Input: the ob Object: OBAtom()
##Output: none
##Return: a bool; is/is not a charged nitrogen
##Calls: none
##Called By: ProTyper(); LigTyper()
def IsChargedN(a):
    atype=a.GetType()
    ChargedN= False
    if not a.IsNitrogen():
        return False
    if atype=='N3':
        num_het=0
        for an in ob.OBAtomAtomIter(a):
            antype=an.GetType()
            if antype!='C3' and (not an.IsHydrogen()):
                num_het+=1
        if num_het==0:
            ChargedN=True
    elif atype=='Ng+':
        if a.ImplicitHydrogenCount()==2 or NumConnectH(a)==2:
            ChargedN=True
    elif a.KBOSum()==4:
        ChargedN=True
    return ChargedN

##-----
##Functions: decide an oxygen atom is a charged oxygen or not
##Input: ob Objece; OBAtom()
##Output: none
##Return: a bool; is/is not a charged oxygen
##Calls: none
##Called By: ProTyper(); LigTyper()
def IsChargedO(a):
    atype=a.GetType()
    ChargedO=False
    if atype in ["O-","O.co2","OCO2"]:
        ChargedO=True
    return ChargedO

def LigAtomTyper(a):
    pmf_type=''
    num_connected_het=a.GetHeteroValence()
    atype=a.GetType()
    if a.IsCarbon():#CF,CP,cF,cP,C3,CW,CO,CN,C0
        if a.IsAromatic() :
            if num_connected_het>0:
                pmf_type='AP'
            else:
                pmf_type="AF"
        else:
            if a.GetHyb()==2:
                if num_connected_het>0:
                    pmf_type="CW"
                else:
                    pmf_type='C3'
            elif a.GetHyb()==1:
                pmf_type='C0'
            else:
                if num_connected_het>0:
                    pmf_type='CP'
                else:
                    pmf_type='CF'
        for an in ob.OBAtomAtomIter(a):
            if IsChargedN(an):
                pmf_type='CN'
                break
        for an in ob.OBAtomAtomIter(a):
            if IsChargedO(an):
                pmf_type='CO'
                break
    elif a.IsOxygen():#OC,OA,OD,OS,OE
        pmf_type='OA'
        if a.IsInRing():
            pmf_type='OE'
        else:
            num_connected_carbon=0
            for an in ob.OBAtomAtomIter(a):
                if an.IsCarbon():
                    num_connected_carbon+=1
            if num_connected_carbon==2:
                pmf_type='OE'
            for an in ob.OBAtomAtomIter(a):
                if not an.IsCarbon() and not an.IsHydrogen():
                    pmf_type='OS'
                    break
        if a.IsHbondDonor() or a.ImplicitHydrogenCount()>0 or NumConnectH(a)>0:
            pmf_type='OD'
        if IsChargedO(a):
            pmf_type='OC'
    elif a.IsNitrogen():#NC,ND,NA,NP,NS,NR,N0
        pmf_type='NA'
        num_connected_h=a.ImplicitHydrogenCount()
        num_connected_c=0
        num_connected_het=0
        for an in ob.OBAtomAtomIter(a):
            if an.IsHydrogen():
                num_connected_h+=1
            elif an.IsCarbon():
                num_connected_c+=1
            else:
                num_connected_het+=1
        if num_connected_c==1 and a.GetHyb()==1:
            pmf_type='N0'
        if num_connected_het>0:
            pmf_type='NS'
        if not a.IsInRing():
            if a.IsHbondAcceptor() or (num_connected_h==0):
                pmf_type='NA'
            if a.IsHbondDonor() or num_connected_h>0:
                pmf_type='ND'
        else:
            if a.GetHyb()==2:
                pmf_type='NR'
                if num_connected_c>=2 and num_connected_h==0 and (not a.IsAromatic()):
                    pmf_type='NP'
        if IsChargedN(a):
            pmf_type='NC'
    elif a.IsSulfur():#SA,SD,SO
        pmf_type='SA'
        num_connected_o=0
        num_connected_h=0
        for an in ob.OBAtomAtomIter(a):
            if an.IsOxygen():
                num_connected_o+=1
            if an.IsHydrogen():
                num_connected_h+=1
        if num_connected_o>0:
            pmf_type='SO'
        if a.IsHbondDonor() or a.ImplicitHydrogenCount()>0 or num_connected_h>0:
            pmf_type='SD'
    elif a.IsHydrogen():
        pmf_type='HL'
    elif a.IsPhosphorus():
        pmf_type='P'
    elif a.GetAtomicNum()==9:
        pmf_type='F'
    elif a.GetAtomicNum()==53:
        pmf_type='I'
    elif a.GetAtomicNum()==35:
        pmf_type='Br'
    elif a.GetAtomicNum()==17:
        pmf_type='CL'
    return pmf_type

def ProAtomTyper(a):
    idx=a.GetIdx()-1
    res=a.GetResidue().GetName()
    atype=a.GetType()
    nidx=a.GetAtomicNum()#Zn,Mg,Ca,Fe,Mn,K
    pmf_type=''
    if a.IsCarbon():#CF,CP,cF,cP,CN,CO
        num_connected_het=0
        for an in ob.OBAtomAtomIter(a):
            antype=an.GetType()
            if not an.IsCarbon() and not an.IsHydrogen():
                num_connected_het+=1
        if not a.IsAromatic():
            if num_connected_het==0:
                pmf_type='CF'
            else:
                pmf_type='CP'
        else:
            if num_connected_het==0:
                pmf_type='AF'
            else:
                pmf_type='AP'
        for an in ob.OBAtomAtomIter(a):
            if IsChargedN(an):
                pmf_type='CN'
                break
        for an in ob.OBAtomAtomIter(a):
            if IsChargedO(an):
                pmf_type='CO'
                break
    elif a.IsOxygen():#OC,OA,OD,OW
        pmf_type='OA'
        if a.IsHbondDonor() or a.ImplicitHydrogenCount()>0 or NumConnectH(a)>0:
            pmf_type='OD'
        if IsChargedO(a):
            pmf_type='OC'
        if res=='HOH' or res=="WAT":
            pmf_type='OW'
    elif a.IsNitrogen():#NC,ND,NA
        pmf_type='NA'
        if a.IsHbondDonor() or a.ImplicitHydrogenCount()>0 or NumConnectH(a)>0:
            pmf_type='ND'
        if IsChargedN(a):
            pmf_type='NC'
    elif a.IsSulfur():#SA,SD
        pmf_type='SA'
        if a.IsHbondDonor() or a.ImplicitHydrogenCount()>0 or NumConnectH(a)>0:
            pmf_type='SD'
    elif a.IsHydrogen():
        pmf_type="HH"
    elif nidx in [30,12,20,26,25,19,11,27,28,29,80,38,48]:#Zn,Mg,Ca,Fe,Mg,K,Na,Co,Ni,Cu,Hg,Sr,Cd
        pmf_type='ME'
    return pmf_type
