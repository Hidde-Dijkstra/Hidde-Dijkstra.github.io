import re

def read_file(file_name):
    f = open(file_name, "r")
    file_content = f.read()
    f.close()
    return file_content

def write_file(file_name, file_content):
    f = open(file_name, 'w')
    f.write(file_content)
    f.close()

def empty(item):
    return item==[] or re.sub(r'[ \n]', '', str(item))==''

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

def join_grid_point(grid_point):
    return '   '.join(['{:.9f}'.format(item)[0:9] for item in grid_point])

def join_grid(k_grid, weight=True):
    grid_str = ''
    for grid_point in k_grid:
        grid_str += '  ' + join_grid_point(grid_point) + ' 1'*weight+'\n'
    return grid_str