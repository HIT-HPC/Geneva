import re
import random

with open('./Patents-format.g', 'w') as wf:
    for i in range(6009554):
        wf.write("v {} {}\n".format(i, random.randint(1, 7)))
    with open('./Patents.g', 'r') as rf:
        line = rf.readline()
        num = 0
        while line:
            infos = re.split(r'[\s]', line)
            if infos[0] == 'v':
                line = rf.readline()
                continue
            if infos[0] == 'e':
                num += 1
                wf.write("{} {} {} {}\n".format(infos[0], infos[1], infos[2], infos[3]))
            line = rf.readline()
print('edge num:{}'.format(num))

