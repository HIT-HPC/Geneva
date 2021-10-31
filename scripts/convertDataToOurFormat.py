import sys
import random
import re
import os

filename = sys.argv[1]
edgelablenum = int(sys.argv[2])

# filename = "DD"
# edgelablenum = 5

path = "./datasets/"
hasedgelabelfile = os.path.exists(path + filename + ".link_labels")
if hasedgelabelfile:
    fedgelabel = open(path + filename + ".link_labels", "r")

hasNodeLabelFile = os.path.exists(path + filename + ".node_labels")
if hasNodeLabelFile:
    fnodel = open(path + filename + ".node_labels", "r")

fedge = open(path + filename + ".edges", "r")
fout = open(path + filename + ".g", "w")
line = fnodel.readline()
unvisitedV = {}
vidindex = 1
while line:
    line = line.strip()
    infos = re.split("[, \t]", line)
    #infos = line.split(" ")
    if infos[0][0] == '#':
        line = fedge.readline()
        continue
    if len(infos) > 1:
        vid = int(infos[0]) - 1
        label = infos[1]
        unvisitedV[vid + 1] = 1
        outstr = "v " + str(vid) + " " + label + "\n"
        fout.write(outstr)
    elif len(infos) == 1:
        label = infos[0]
        unvisitedV[vidindex] = 1
        outstr = "v " + str(vidindex) + " " + label + "\n"
        fout.write(outstr)
    line = fnodel.readline()
    vidindex = vidindex + 1

line = fedge.readline()
if hasedgelabelfile:
    edgelabel = fedgelabel.readline()
    edgelabel.strip()
visitededge = {}
linenum = 1
lastvid = 0
while line:
    line = line.strip()
    line = line.replace(" ", "")
    infos = re.split(r'[,\s]', line)
    key1 = infos[0] + " " + infos[1]
    key2 = infos[1] + " " + infos[0]
    # print(infos[0])
    # print(infos[1])
    # print(key1, key2)
    if (key1 in visitededge) or (key2 in visitededge):
        line = fedge.readline()
        linenum = linenum + 1
        continue
    visitededge[key1] = linenum
    visitededge[key2] = linenum
    svid = int(infos[0]) - 1
    evid = int(infos[1]) - 1
    lastvid = evid + 1
    if svid + 1 in unvisitedV:
        del unvisitedV[svid + 1]
    if evid + 1 in unvisitedV:
        del unvisitedV[evid + 1]
    if not hasedgelabelfile:
        edgelabel = str(random.randint(0, edgelablenum - 1))
    outstr = "e " + str(svid) + " " + str(evid) + " " + edgelabel + "\n"
    fout.write(outstr)
    line = fedge.readline()
    if hasedgelabelfile:
        edgelabel = fedgelabel.readline()
        edgelabel.strip()
    linenum = linenum + 1
if len(unvisitedV) > 0:
    print(str(len(unvisitedV)) + " vertices are unvisited by edges")
    keys = list(unvisitedV.keys())
    firstkey = keys[0]
    outstr = "e " + str(lastvid - 1) + " " + str(firstkey - 1) + " 0\n"
    fout.write(outstr)
    for i in range(1, len(keys)):
        key = keys[i]
        outstr = "e " + str(firstkey - 1) + " " + str(key - 1) + " 0\n"
        fout.write(outstr)

fedge.close()
fnodel.close()
fout.close()
