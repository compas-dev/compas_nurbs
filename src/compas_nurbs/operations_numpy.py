import numpy as np


def normalize(array):
    array = np.array(array)
    return (array.T/np.linalg.norm(array, axis=1).reshape(1, -1)).T


def surface_normals(surface, params):
    skl = surface.derivatives_at(params, order=1)
    vector = np.cross(skl[:, 1, 0], skl[:, 0, 1])
    return normalize(vector)
