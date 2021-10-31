import os
import sys
import random
import re
from logger import Logger


def readNodeAndWriteToFile(path, fileName, nodeLabelNum):
    hasNodeLabelFile = os.path.exists(path + fileName + ".node_labels")
    fout = open(path + fileName + ".g", "w")
    unvisitedV = {}
    nodeDict = {}
    if hasNodeLabelFile:
        fnodeLabel = open(path + fileName + ".node_labels", "r")
        line = fnodeLabel.readline()
        vidindex = 1
        nodeNum = 0
        while line:
            line = line.strip()
            infos = re.split(r'[, \s]', line)
            if infos[0][0] == '#':
                line = fnodeLabel.readline()
                continue
            if len(infos) > 1:
                vid = int(infos[0])
                label = infos[1]
                unvisitedV[vid] = 1
                nodeDict.update({vid: 0})
                outstr = "v {} {}\n".format(str(vid), label)
                fout.write(outstr)
            elif len(infos) == 1:
                label = infos[0]
                unvisitedV[vidindex] = 1
                nodeDict.update({vidindex: 0})
                outstr = "v {} {}\n".format(str(vidindex), label)
                fout.write(outstr)
            line = fnodeLabel.readline()
            vidindex = vidindex + 1
            nodeNum = nodeNum + 1
        fnodeLabel.close()
        fout.close()
        print("the number of {}\'s node: {}.".format(fileName, nodeNum))
    else:
        hasEdgeFile = os.path.exists(path + fileName + ".edges")
        if hasEdgeFile:
            fedge = open(path + fileName + ".edges", "r")
        else:
            fedge = open(path + fileName + ".txt", "r")
        line = fedge.readline()
        while line:
            infos = re.split(r'[, \s]', line)
            if infos[0][0] == '#' or infos[0][0] == '%':
                line = fedge.readline()
                continue
            nodeDict.update({int(infos[0]): 0})
            nodeDict.update({int(infos[1]): 0})
            line = fedge.readline()
        nodes = list(nodeDict.keys())
        print('the number of {}\'s node: {}.'.format(fileName, len(nodes)))
        for i in range(len(nodes)):
            unvisitedV[nodes[i]] = 1
            outstr = "v {} {}\n".format(nodes[i], random.randint(1, nodeLabelNum))
            fout.write(outstr)
        fedge.close()
        fout.close()
    return unvisitedV, nodeDict


def readEdgeAndWriteToFile(path, fileName, edgeLableNum, unvisitedV, nodeDict):
    hasEdgeLabelFile = os.path.exists(path + fileName + ".link_labels")
    if hasEdgeLabelFile:
        fedgeLabel = open(path + fileName + ".link_labels", "r")
    hasEdgeFile = os.path.exists(path + fileName + ".edges")
    if hasEdgeFile:
        fedge = open(path + fileName + ".edges", "r")
    else:
        fedge = open(path + fileName + ".txt")
    fout = open(path + fileName + ".g", "a")
    line = fedge.readline()
    if hasEdgeLabelFile:
        edgeLabel = fedgeLabel.readline()
        edgeLabel.strip()
    visitedEdge = {}
    lineNum = 1
    lastVid = 0
    edgeNum = 0
    while line:
        line = line.replace(" ", "")
        infos = re.split(r'[, \s]', line)
        if infos[0][0] == '#' or infos[0][0] == '%':
            line = fedge.readline()
            continue
        key1 = "{} {}".format(infos[0], infos[1])
        key2 = "{} {}".format(infos[1], infos[0])
        if (key1 in visitedEdge) or (key2 in visitedEdge):
            line = fedge.readline()
            lineNum = lineNum + 1
            continue
        visitedEdge[key1] = lineNum
        visitedEdge[key2] = lineNum
        sVid = int(infos[0])
        eVid = int(infos[1])
        lastvid = eVid
        if sVid in unvisitedV:
            del unvisitedV[sVid]
        if eVid in unvisitedV:
            del unvisitedV[eVid]
        if not hasEdgeLabelFile:
            edgeLabel = str(random.randint(1, edgeLableNum))
        outstr = "e {} {} {}\n".format(sVid, eVid, edgeLabel)
        fout.write(outstr)
        nodeDict[sVid] = nodeDict[sVid] + 1
        nodeDict[eVid] = nodeDict[eVid] + 1
        edgeNum = edgeNum + 1
        line = fedge.readline()
        if hasEdgeLabelFile:
            edgeLabel = fedgeLabel.readline()
            edgeLabel.strip()
        lineNum = lineNum + 1
    if len(unvisitedV) > 0:
        print("{} vertices are unvisited by edges.".format(len(unvisitedV)))
        keys = list(unvisitedV.keys())
        firstKey = keys[0]
        outstr = "e {} {} 0\n".format(lastVid, firstKey)
        fout.write(outstr)
        for i in range(1, len(keys)):
            outstr = "e {} {} 0\n".format(firstKey, keys[i])
            fout.write(outstr)
    fedge.close()
    fout.close()
    maxDegree = max(nodeDict.values())
    sumDegree = sum(nodeDict.values())
    averageDegree = sumDegree / len(nodeDict.keys())
    print("the number of {}\'s edge: {}.".format(fileName, edgeNum))
    print("the max degree of {} is {}, the average degree of {} is {:.2f}.".
          format(fileName, maxDegree, fileName, averageDegree))


if __name__ == '__main__':
    sys.stdout = Logger(sys.stdout)
    sys.stderr = Logger(sys.stderr)
    path = "./datasets/"
    fileName = sys.argv[1]
    nodeLabelNum = int(sys.argv[2])
    edgeLableNum = int(sys.argv[3])
    unvisitedV, nodeDict = readNodeAndWriteToFile(path, fileName, nodeLabelNum)
    readEdgeAndWriteToFile(path, fileName, edgeLableNum, unvisitedV, nodeDict)

