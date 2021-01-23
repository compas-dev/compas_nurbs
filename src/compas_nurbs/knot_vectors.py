import math
from compas.geometry import Bezier
from compas.geometry import distance_point_point


class CurveKnotStyle(object):
    """
    """
    Uniform = 0
    Chord = 1
    ChordSquareRoot = 2

    def __init__(self, style, periodic=False):
        self.style = style
        self.periodic = periodic

    # TODO add methods ?


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


def bezier_spaced_parameters(num_points):
    """Computes bezier spaced parameters for interpolation.

    Parameters
    ----------
    num_points : int
        Number of points to compute parameters for.

    Returns
    -------
    list of float
        Parameters within the range of [0, 1].

    Notes
    -----
    If degree >= 5 this provides a more stable knotvector for 'uniform' knot style.
    """
    curve = Bezier([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 1.0, 0.0], [1., 1., 0]])
    return [curve.point(t).y for t in equally_spaced_parameters(num_points)]


def chord_lengths(points):
    """Returns the chord lengths.

    Parameters
    ----------
    points : list of point
        Points on the curve.

    Returns
    -------
    list of float
        The n - 1 lengths between the n points.

    Examples
    --------
    >>>
    """
    return [distance_point_point(a, b) for a, b in zip(points[:-1], points[1:])]


def chord_spaced_parameters(points):
    """Computes parameters based on chord length for interpolation.

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

    Examples
    --------
    >>>
    """
    cds = chord_lengths(points)
    d = sum(cds)
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

    Examples
    --------
    >>>
    """
    cds = [math.sqrt(v) for v in chord_lengths(points)]
    d = sum(cds)
    uk = [0.] + cds + [1.]
    return [sum(uk[0:i + 1]) / d for i in range(len(points))]


def knot_vector_uniform(num_points, degree, periodic=False):
    """Computes a uniform knot vector.

    Parameters
    ----------
    num_points : int
        Number of points to compute parameters for.
    degree : int
        The degree of the curve.

    Returns
    -------
    list of float
        The knot vector in the domain of [0, 1].

    Notes
    -----
    Is the same to which Rhino refers to as CurveKnotStyle.Uniform

    """
    kv = [0.0 for _ in range(degree + 1)]
    step = 1. / (num_points - degree)
    for i in range(1, num_points - degree):
        kv.append(step * i)
    kv += [1.0 for _ in range(degree + 1)]
    return kv


def knot_vector_from_params(degree, params, periodic=False):
    """Computes a knot vector from parameters using the averaging method.

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

    Notes
    -----
    Is the same as geomdl.fitting.compute_knot_vector
    """
    kv = [0.0 for _ in range(degree + 1)]
    for j in range(1, len(params) - degree):
        v = (1.0 / degree) * sum([params[j] for j in range(j, j + degree)])
        kv.append(v)
    kv += [1.0 for _ in range(degree + 1)]
    return kv


def knot_vector_and_params(points, degree, knot_style, extended=False, periodic=False):
    """Returns the knot vector for curve interpolation.

    Parameters
    ----------
    points : list of point
        Points on the curve.
    degree : int
        The degree of the curve
    knot_style : int
        The knot style: 0 for 'Uniform', 1 for 'Chord', 2 for 'ChordSquareRoot'
    extended : bool
        `True` if the knot vector and params should be calculated for interpolation
        with end derivatives. Defaults to `False`.
    periodic : bool
        `True` if knot vector should be computed for a `closed` curve. Defaults to
        `False`.

    Returns
    -------
    tuple : (knot_vector, parameters)
        The knot vector and the parameters for the interpolation.
    """

    if knot_style == CurveKnotStyle.Uniform:
        uk = equally_spaced_parameters(len(points))
    elif knot_style == CurveKnotStyle.Chord:
        uk = chord_spaced_parameters(points)
    elif knot_style == CurveKnotStyle.ChordSquareRoot:
        uk = chord_sqrt_spaced_parameters(points)
    else:
        raise ValueError("CurveKnotStyle %d is unknown" % knot_style)

    if extended:  # extend parameters for end derivatives estimation
        uk = [uk[0], 0.] + uk[1:-1] + [1., uk[-1]]

    # if knot_style == CurveKnotStyle.Uniform:
    #    kv = knot_vector_uniform(len(uk), degree)
    # else:
    kv = knot_vector_from_params(degree, uk)

    return kv, uk


def normalize_knot_vector(knot_vector):
    """Returns a normalized knot vector within the [0, 1] domain.

    Parameters
    ----------
    list of float
        A knot vector

    Returns
    -------
    list of float
        The normalized knot vector.
    """
    knot_vector = [v - knot_vector[0] for v in knot_vector]
    return [v / float(knot_vector[-1]) for v in knot_vector]


def check_knot_vector(knot_vector, num_points, degree):
    """Checks the validity of the knot vector.

    Parameters
    ----------
    list of float
        The knot vector to check.
    num_points : int
        Number of points to compute parameters for.
    degree : int
        The degree

    Returns
    -------
    bool
        True if the knot vector is valid, False otherwise.
    """
    if len(knot_vector) != degree + num_points + 1:  # m = p + n + 1
        return False

    # Check ascending order
    prev_knot = knot_vector[0]
    for knot in knot_vector:
        if prev_knot > knot:
            return False
        prev_knot = knot
    return True
