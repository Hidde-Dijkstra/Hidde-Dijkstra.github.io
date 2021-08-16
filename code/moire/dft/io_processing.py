import re
import subprocess
from .decorators import time

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
        return [float(re.sub("[^0-9.-]", "", x)) 
                for x in string.split(char)]
    
def delim(string, lim1, lim2, lim3=None):
    if lim3 == None:
        return string.split(lim1)[1].split(lim2)[0]
    else:
        return string.split(lim1)[1].split(lim2)[0].split(lim3)[0]

def gen_lst(lst, str, func=lambda x: x, ignore_first=False):
    return [func(item) for i, item in enumerate(lst.split(str))
            if (not empty(item)) and (i!=0 or not ignore_first)]

def join_grid_point(grid_point):
    return '   '.join(['{:.9f}'.format(item)[0:9] for item in grid_point])

def join_grid(k_grid, weight=True):
    return '\n'.join([f"  {join_grid_point(grid_point)}{' 1'*weight}"
                      for grid_point in k_grid])
 
@time   
def run(command=''): 
    print(command)
    print(subprocess.run(command, 
                         capture_output=True, 
                         text=True,
                         shell=True).stdout)