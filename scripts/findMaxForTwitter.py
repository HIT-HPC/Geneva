import re
import os

with open("./twitter/twitter-uniq.txt", "r") as f:
    nodeSet = set()
    line = f.readline()
    while line:
        infos = re.split(r'[\s]', line)
        nodeSet.add(int(infos[0]))
        nodeSet.add(int(infos[1]))
        line = f.readline()
sortSet = sorted(nodeSet)
with open("./twitter.g", "w") as wf:
    for i in sortSet:
        wf.write("v {}\n".format(i))


