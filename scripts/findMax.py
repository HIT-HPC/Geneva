import re
import os
import sys

file = sys.argv[1]

with open("./{}-edges.g".format(file), "r") as f:
    maxLabel = 0
    minLabel = 999999999
    num = 0
    line = f.readline()
    while line:
        infos = re.split(r'[, \s]', line)
        infos = [item for item in filter(lambda x:x != '', infos)]
        if infos[0] == '#' or infos[0] == '%':
            line = f.readline()
            continue
        num += 1
        maxLabel = max(int(infos[0]), int(infos[1]), maxLabel)
        minLabel = min(int(infos[0]), int(infos[1]), minLabel)
        line = f.readline()
    print("max: {}, min: {}, num: {}".format(maxLabel, minLabel, num))


