 #!/usr/bin/env python
#-------------------------------------------------------------------------------
# Name: Format.py
# Purpose: Generate data in Libsvm format
# Author:  BearShare
# Created:     21-05-2013
# Licence:     <your licence>
# Last Modified: 
#-------------------------------------------------------------------------------

import os,string

def GetList():
    w=[]
    f=open("list")
    for line in f.readlines():
        w.append(line.rstrip().lstrip())
    return w
def main():
    fo=open("LibSvm_GlideCfingerv.txt",'wb')
    file_list=GetList()
    for each in file_list:
        if os.path.isfile("Cfingerv/%s_GlideCfingerv.txt"%each):
            print "%s_GlideCfingerv.txt"%each
            Aline=''
            f1=open("Cfingerv/%s_GlideCfingerv.txt"%each)
            f2=open("ic50/%s.txt"%each)
            Pic50=float(f2.readline().rstrip().lstrip())
            Aline+="%.2f "%Pic50
            lines=f1.readlines()
            for i in range(len(lines)):
                aline=lines[i].rstrip().lstrip().split(" ")
                #print aline
                value=float(aline[-1])
                if value == 0.00:
                    continue
                else:
                    Aline+="%d:%.2f "%(i+1,value)
            Aline+="\n"
            fo.write(Aline)
        else:
            continue
    fo.close()

if __name__=="__main__":
    main()
