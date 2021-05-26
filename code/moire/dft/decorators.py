import os
import numpy as np
from datetime import datetime


def change_directory(directory_name):
    def decorator(func):
        def wrapper(self, *args, **kwargs):
            current_dir = os.getcwd()
            new_dir = self.__dict__[directory_name]
            if not os.path.isdir(new_dir):
                os.mkdir(new_dir)
            os.chdir(new_dir)
            output = func(self, *args, **kwargs)
            os.chdir(current_dir)
            return output
        return wrapper
    return decorator

def time(func):
    def wrapper(*args, **kwargs):
        startTime = datetime.now()
        output = func(*args, **kwargs)
        print(datetime.now()-startTime)
        return output
    return wrapper

def save(func):
    def wrapper(self, *args, **kwargs):
        data_dic = func(self, *args, **kwargs)
        for file_name, array in data_dic.items():
            file_loc = self.prefix+'/'+file_name
            folders = file_loc.split('/')
            for folder in ['/'.join(folders[:i]) for i in range(1, len(folders))]:
                if not os.path.isdir(folder):
                    os.mkdir(folder)           
            np.save(file_loc, np.array(array))
    return wrapper

def check(*check_args):
    def decorator(func):
        def wrapper(self, *args, **kwargs):
            bools = [arg for arg in check_args if not self.__dict__[arg]]
            if len(bools) == 0:
                return func(self, *args, **kwargs)
            else:
                raise Exception(' '.join(bools) + ' not prepared')
        return wrapper
    return decorator
