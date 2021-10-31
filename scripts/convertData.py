import re
import sys

fileName = sys.argv[1]

with open('./'+fileName+'.g', 'r') as rf:
    with open('./'+fileName+'-format.g', 'w') as wf:
        line = rf.readline()
        while line:
            infos = re.split(r'[\s]', line)
            if infos[0] == 'v':
                wf.write('v {} {}\n'.format(int(infos[1]) + 1, infos[2]))
            if infos[0] == 'e':
                wf.write('e {} {} {}\n'.format(int(infos[1]) +1, int(infos[2]) + 1, infos[3]))
            line = rf.readline()
