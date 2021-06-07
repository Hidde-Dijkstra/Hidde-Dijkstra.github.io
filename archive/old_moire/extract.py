from .decorators import change_directory, save
import numpy as np
import re

@change_directory('data_directory')
@save
@change_directory('work_directory')
def extract_relax(self):
    f = read_file(self.prefix+'relax.out')
    position_data = f.split('ATOMIC_POSITIONS (alat)\n')[1:]
    position_data = [item.split('\nEnd final')[0] for item in position_data]
    data_dic = {}
    for i, part in position_data:
        for row in [gen_lst(row, ' ') for row in part.split('\n')]:
            dic_key = 'iter_' + str(i) + '_' + row[0]
            if dic_key in data_dic:
                data_dic[dic_key].append(row[1:])
            else:
                data_dic[dic_key] = [row[1:]]
    return data_dic

def scrub_str(string, char=None):
    if char == None:
        return float(re.sub("[^0-9.-]", "", string))
    else:
        return [float(re.sub("[^0-9.-]", "", x)) for x in string.split(char)]

def gen_lst(lst, str, func=lambda x: x, ignore_first=False):
    new_lst = []
    for i, item in enumerate(lst.split(str)):
        if (not empty(item)) and (i!=0 or not ignore_first):
            new_lst.append(func(item))
    return new_lst

def read_file(file_name):
    f = open(file_name, "r")
    file_content = f.read()
    f.close()
    return file_content

def empty(item):
    return item==[] or re.sub(r'[ \n]', '', str(item))==''