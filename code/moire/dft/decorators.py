import os
import numpy as np
from datetime import datetime


def change_directory(directory_name):
    def decorator(func):
        def wrapper(self, *args, **kwargs):
            current_dir = os.getcwd()
            os.chdir(self.__dict__[directory_name])
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
        for file_name, np_array in data_dic.items():
            np.save(self.prefix+'/'+file_name, np_array)
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
