import sys
import random
import re

def constructDict(filename):
    emb=open(filename,"r")
    line = emb.readline()
    tmpdict={}
    while line:
        line = line.strip()
        if len(line)==0:
            line = emb.readline()
            continue
        if line in tmpdict:
            tmpdict[line] = tmpdict[line]+1
        else:
            tmpdict[line] = 1
        line = emb.readline()
    
    emb.close()
    return tmpdict




gpufilename = sys.argv[1]
cpufilename = sys.argv[2]
fout = open("compareResults.txt","w")
gpudict = constructDict(gpufilename)
cpudict = constructDict(cpufilename)
bothlist=[]
onlygpulistzero=[]
onlygpulistwrong=[]
onlycpulist=[]
morelist=[]
for key in gpudict:
    if gpudict[key]>1:
        morelist.append(key)
    else:
        if key in cpudict:
            bothlist.append(key)
            cpudict.pop(key)
        else:
            eles = key.split(" ")
            found=0
            for vid in eles:
                if vid=="0":
                    found=1
            if found==0:
                onlygpulistwrong.append(key)
            else:
                onlygpulistzero.append(key)

for key in cpudict:
    onlycpulist.append(key)

fout.write("###########in both cpu and gpu ("+str(len(bothlist))+") ###########\n")
for item in bothlist:
    fout.write(item+"\n")
fout.write("###########only in gpu zero ("+str(len(onlygpulistzero))+") ###########\n")
for item in onlygpulistzero:
    fout.write(item+"\n")
fout.write("###########only in gpu wrong ("+str(len(onlygpulistwrong))+") ###########\n")
for item in onlygpulistwrong:
    fout.write(item+"\n")
fout.write("###########only in cpu ("+str(len(onlycpulist))+") ###########\n")
for item in onlycpulist:
    fout.write(item+"\n")
fout.write("###########more than 1 ("+str(len(morelist))+") ###########\n")
for item in morelist:
    fout.write(item+"\n")
fout.close()