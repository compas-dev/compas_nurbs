"""
Checkout
https://github.com/johntfoster/bspline/blob/master/bspline/splinelab.py

 knotstyle(opt):
0 Uniform knots.  Parameter spacing between consecutive knots is 1.0.
1 Chord length spacing.  Requires degree = 3 with arrCV1 and arrCVn1 specified.
2 Sqrt (chord length).  Requires degree = 3 with arrCV1 and arrCVn1 specified.
3 Periodic with uniform spacing.
4 Periodic with chord length spacing.  Requires an odd degree value.
5 Periodic with sqrt (chord length) spacing.  Requires an odd degree value.

"""
import math
from compas.geometry import Bezier
from compas.geometry import distance_point_point

def equally_spaced_parameters(num_points):
    """Computes equally spaced parameters for interpolation.

    Please refer to the equation 9.3 on The NURBS Book (2nd Edition), pp.364.

    Parameters
    ----------
    num_points : int
        Number of points to compute parameters for.

    Returns
    -------
    list of float
        Parameters within the range of [0, 1].
    """
    n = num_points - 1
    return [i / float(n) for i in range(n + 1)]

def chord_lengths(points):
    """Computes the 
    """
    #return np.linalg.norm(points[:-1] - points[1:], axis=1)
    return [distance_point_point(a, b) for a, b in zip(points[:-1], points[1:])]

def chord_spaced_parameters(points):
    """Computes parameters based on the chord length for interpolation.

    Please refer to the Equations 9.4 and 9.5 for chord length parametrization,
    and Equation 9.6 for centripetal method on The NURBS Book (2nd Edition),
    pp.364-365.

    Parameters
    ----------
    points : list of point
        Points on the curve.

    Returns
    -------
    list of float
        Parameters within the range of [0, 1].

    Note
    ----
    This spacing is actually the one to which Rhino refers to as `CurveKnotStyle.Uniform`.

    Examples
    --------
    >>>
    """
    cds = chord_lengths(points)
    d = sum(cds)  # total chord length
    #uk = np.append(np.append([0], cds), [1.])
    uk = [0.] + cds + [1.]
    return [sum(uk[0:i + 1]) / d for i in range(len(points))]


def chord_sqrt_spaced_parameters(points):
    """Computes parameters based on the chord length for interpolation.

    Please refer to the Equations 9.4 and 9.5 for chord length parametrization,
    and Equation 9.6 for centripetal method on The NURBS Book (2nd Edition),
    pp.364-365.

    Parameters
    ----------
    points : list of point
        Points on the curve.

    Returns
    -------
    list of float
        Parameters within the range of [0, 1].

    Note
    ----
    This spacing is actually the one to which Rhino refers to as `CurveKnotStyle.Uniform`.

    Examples
    --------
    >>>
    """
    #cds = np.sqrt(chord_lengths(points))
    cds = [math.sqrt(v) for v in chord_lengths(points)]
    d = sum(cds)
    uk = [0.] + cds + [1.]
    return [sum(uk[0:i + 1]) / d for i in range(len(points))]





def bezier_spaced_parameters(num_points):
    """
    """
    curve = Bezier([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 1.0, 0.0], [1., 1., 0]])
    return [curve.point(t).y for t in equally_spaced_parameters(num_points)]


def knot_vector_uniform(num_points, degree, periodic=False):
    """
    """ 
    kv = [0.0 for _ in range(degree + 1)]
    step = 1. / (num_points - degree)
    for i in range(1, num_points - degree):
        kv.append(step * i)
    kv += [1.0 for _ in range(degree + 1)]
    return kv


def knot_vector_from_params(degree, params, periodic=False):
    """ Computes a knot vector from parameters using the averaging method.

    Please refer to the Equation 9.8 on The NURBS Book (2nd Edition), pp.365 for
    details.

    Parameters
    ----------
    degree : int
        The degree of the curve
    params : list of float
        Parameters on the curve in the range of [0, 1].

    Returns
    -------
    list of float
        The knot vector.

    == from geomdl.fitting.compute_knot_vector
    """
    kv = [0.0 for _ in range(degree + 1)]
    for j in range(1, len(params) - degree):
        v = (1.0 / degree) * sum([params[j] for j in range(j, j + degree)])
        kv.append(v)
    kv += [1.0 for _ in range(degree + 1)]
    return kv


def knot_vector_and_params(points, degree, knot_style, extended=False, periodic=False):
    """
    """

    if knot_style == 'uniform':
        uk = equally_spaced_parameters(len(points)) if degree <=3 else bezier_spaced_parameters(len(points))
    elif knot_style == 'chord':
        uk = chord_spaced_parameters(points)
    elif knot_style == 'chord_sqrt':
        uk = chord_sqrt_spaced_parameters(points)
    else:
        raise NotImplementedError

    if extended: # extend parameters for end derivatives estimation
        uk = [uk[0], 0.] + uk[1:-1] + [1., uk[-1]]

    if knot_style == 'uniform':
        kv = knot_vector_uniform(len(uk), degree)
    else:
        kv = knot_vector_from_params(degree, uk)

    return kv, uk


def normalize_knot_vector(knot_vector):
    """Normalizes the input knot vector to [0, 1] domain.
    """
    k = float(knot_vector[-1])
    return [v/k for v in knot_vector]
