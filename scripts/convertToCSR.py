'''
Author: your name
Date: 2021-08-28 10:33:18
LastEditTime: 2021-08-28 12:20:04
LastEditors: Please set LastEditors
Description: In User Settings Edit
FilePath: \convertToCSR\convertToCSR.py
'''

import re
import sys
import struct

file = sys.argv[1]
row_num = int(sys.argv[2])
value_num = int(sys.argv[3])

# file = "test"
# row_num = 4
# value_num = 9

values = [0 for i in range(value_num)]
col_num = [0 for i in range(value_num)]
row_offset = [0 for i in range(row_num + 1)]
value_offset = 0
edge_sum = 0
current_start_sum = 0
last_row_num = 0
min_col_num = 9999999999

with open('./{}-edges.g'.format(file), 'r') as rf:
    line = rf.readline()
    while line:
        infos = re.split(r'[\s]', line)
        infos = [item for item in filter(lambda x:x != '', infos)]
        s_num = int(infos[0])
        e_num = int(infos[1])
        label = int(infos[2])
        if e_num > s_num:
            s_num, e_num = e_num, s_num
        if last_row_num == s_num:
            min_col_num = min(min_col_num, e_num)
        elif last_row_num < s_num:
            row_offset[last_row_num] = edge_sum
            min_col_num = 9999999999
            min_col_num = min(min_col_num, e_num)
            last_row_num = s_num
            edge_sum += current_start_sum
            current_start_sum = 0
        values[value_offset] = label
        col_num[value_offset] = e_num
        value_offset += 1
        current_start_sum += 1
        line = rf.readline()
    row_offset[last_row_num] = edge_sum
    row_offset[row_num] = value_num

# print(row_offset)
# print(col_num)
# print(values)
print("row_offset size:{}".format(sys.getsizeof(row_offset)))
print("col_num size:{}".format(sys.getsizeof(col_num)))
print("values size:{}".format(sys.getsizeof(values)))

with open('./{}.csr'.format(file), 'wb') as wf:
    for i in row_offset:
        s = struct.pack('I', i)
        wf.write(s)
    for i in col_num:
        s = struct.pack('I', i)
        wf.write(s)
    for i in values:
        s = struct.pack('I', i)
        wf.write(s)
print("done")

