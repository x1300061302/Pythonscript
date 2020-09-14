from const import *
import h5py
def get_f(name):
    f = h5py.File(name)
    return f;

