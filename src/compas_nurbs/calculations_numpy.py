import numpy as np


def normalize(array):
    array = np.array(array)
    return (array.T/np.linalg.norm(array, axis=1).reshape(1, -1)).T
