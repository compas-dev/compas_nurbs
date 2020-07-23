"""
Checkout
https://github.com/johntfoster/bspline/blob/master/bspline/splinelab.py
"""


def knot_vector_uniform(degree, number_of_control_points, open=True, periodic=False):
    """
    """
    if number_of_control_points - degree < 2:
        raise ValueError("Degree must be smaller or equal to %d." % (number_of_control_points - 2))

    kv = [0.0 for _ in range(degree + 1)]
    step = 1. / (number_of_control_points - degree)
    for i in range(1, number_of_control_points - degree):
        kv.append(step * i)
    kv += [1.0 for _ in range(degree + 1)]
    return kv


def knot_vector_chord(degree, periodic=False):
    """
    """
    raise NotImplementedError


def knot_vector_chord_square_root(degree, periodic=False):
    """
    """
    raise NotImplementedError
