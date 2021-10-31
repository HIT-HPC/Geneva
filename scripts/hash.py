import sys
import hashlib

def gen_md5(data):
    md5 = hashlib.md5()
    md5.update(data.encode('utf-8'))
    return md5.hexdigest()
    
def big_file_remove_same(input_file, output_file):
    finger_print_set = set()
    with open(input_file, 'r', encoding='utf-8') as f, open(output_file, 'w', encoding='utf-8') as ff,open('dumple.txt', 'w', encoding='utf-8') as fff:
        for line in f:
            finger_print = gen_md5(line)
            if finger_print not in finger_print_set:
                finger_print_set.add(finger_print)
                ff.write(line)
            else:
                fff.write(line)


if __name__ == "__main__":
    big_file_remove_same(sys.argv[1],sys.argv[2])

