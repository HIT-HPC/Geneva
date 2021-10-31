'''
Author: your name
Date: 2021-08-28 10:21:44
LastEditTime: 2021-08-28 10:30:28
LastEditors: Please set LastEditors
Description: In User Settings Edit
FilePath: \convertToCSR\extractEdges.py
'''
import re
import sys

file = sys.argv[1]

with open('./{}.g'.format(file), 'r') as rf:
    with open('./{}-edges.g'.format(file), 'w') as wf:
        line = rf.readline()
        vnum = 0
        enum = 0
        while line:
            infos = re.split(r'[\s]', line)
            if infos[0] == 'v':
                vnum += 1
                line = rf.readline()
                continue
            enum += 1
            wf.write("{} {} {}\n".format(infos[1], infos[2], infos[3]))
            line = rf.readline()
            
