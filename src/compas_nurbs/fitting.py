import numpy as np
from scipy.linalg import lu_factor, lu_solve

from compas.geometry import allclose  # noqa: F401

from compas_nurbs.helpers import find_spans
from compas_nurbs.helpers import basis_functions
from compas_nurbs.knot_vectors import knot_vector_and_params
from compas_nurbs.knot_vectors import CurveKnotStyle


# TODO: estimate derivatives for degree==3
# https://link.springer.com/content/pdf/10.1007/s003660050038.pdf

def global_curve_interpolation(points, degree, knot_style=0, periodic=False):
    """Global curve interpolation through points.

    Please refer to Algorithm A9.1 on The NURBS Book (2nd Edition), pp.369-370
    for details.

    Parameters
    ----------
    points : list of point
        The list of points on the curve we are looking for.
    degree : int
        The degree of the output parametric curve.
    knotstyle : int, optional
        The knot style, either 0, 1, or 2 [uniform, chord, or chord_square_root].
        Defaults to 0, uniform.

    Returns
    -------
    tuple (control_points, knot_vector)

    Examples
    --------
    >>> degree = 3
    >>> points = [(0., 0., 0.), (3., 4., 0.), (-1., 4., 0.), (-4., 0., 0.), (-4., -3., 0.)]
    >>> control_points, knot_vector = global_curve_interpolation(points, degree)
    >>> solution = [(0, 0, 0), (6.44, 3.72, 0.), (-2.67, 7.5, 0), (-5.11, -2.72, 0), (-4, -3, 0)]
    >>> allclose(control_points, solution, tol=0.005)
    True
    """
    points = np.array(points)
    kv, uk = knot_vector_and_params(points, degree, knot_style, extended=False)
    M = coefficient_matrix(degree, kv, uk)
    lu, piv = lu_factor(M)
    return lu_solve((lu, piv), points), kv


def global_curve_interpolation_with_end_derivatives(points,
                                                    degree,
                                                    start_derivative,
                                                    end_derivative,
                                                    knot_style=0,
                                                    periodic=False):
    """Global curve interpolation through points with end derivatives specified.

    Please refer to The NURBS Book (2nd Edition), pp. 370-371 for details.

    Parameters
    ----------
    points : list of point
        The list of points on the curve we are looking for.
    degree : int
        The degree of the output parametric curve.
    start_derivative : vector
        The start derivative of the curve
    end_derivative : vector
        The end derivative of the curve
    knotstyle : int, optional
        The knot style, either 0, 1, or 2 [uniform, chord, or chord_square_root].
        Defaults to 0, uniform.

    Returns
    -------
    tuple (control_points, knot_vector)

    Examples
    --------
    >>> degree = 3
    >>> points = [(0., 0., 0.), (3., 4., 0.), (-1., 4., 0.), (-4., 0., 0.), (-4., -3., 0.)]
    >>> start_derivative = [17.75, 10.79, 0.0]
    >>> end_derivative = [-0.71, -12.62, 0.0]
    >>> result = global_curve_interpolation_with_end_derivatives(points, degree, start_derivative, end_derivative)
    >>> solution = [(0, 0, 0), (1.48, 0.9, 0), (4.95, 5.08, 0), (-1.56, 4.87, 0), (-4.72, -0.56, 0), (-3.94, -1.95, 0), (-4, -3, 0)]
    >>> allclose(result[0], solution, tol=0.005)
    True
    """

    points = np.array(points)
    kv, uk = knot_vector_and_params(points, degree, knot_style, extended=True)
    M = coefficient_matrix(degree, kv, uk)

    M[1][0] = -1.
    M[1][1] = 1.
    M[-2][-2] = -1.
    M[-2][-1] = 1.

    v0 = np.array(start_derivative) * kv[degree + 1] / degree
    vn = np.array(end_derivative) * (1 - kv[len(kv) - 1 - degree - 1]) / degree
    C = points[:]
    C = np.insert(C, 1, v0, axis=0)
    C = np.insert(C, -1, vn, axis=0)

    lu, piv = lu_factor(M)
    return lu_solve((lu, piv), C), kv


def coefficient_matrix(degree, knot_vector, uk):
    """Returns the coefficient matrix for global interpolation.

    Please refer to The NURBS Book (2nd Edition), pp. 370-371 for details.

    Parameters
    ----------
    degree : int
        The degree of the curve.
    knot_vector : list of float
        The knot vector of the curve.
    uk : list of float
        parameters

    Returns
    -------
    coefficient matrix, list
    """
    # TODO: use this for evaluators?
    num_points = len(uk)
    spans = find_spans(knot_vector, num_points, uk)
    bases = basis_functions(degree, knot_vector, spans, uk)
    M = [[0.0 for _ in range(num_points)] for _ in range(num_points)]
    for i, (span, basis) in enumerate(zip(spans, bases)):
        M[i][span - degree:span + 1] = basis[:degree + 1]
    return np.array(M)


def interpolate_curve(points, degree, knot_style=0, start_derivative=None, end_derivative=None, periodic=False):
    """Interpolate curve by the specified parameters.

    Parameters
    ----------
    points : list of point
        The list of points on the curve we are looking for.
    degree : int
        The degree of the output parametric curve.
    start_derivative : vector, optional
        The start derivative of the curve. Defaults to ``None``.
    end_derivative : vector
        The end derivative of the curve. Defaults to ``None``.
    knotstyle : int, optional
        The knot style, either 0, 1, or 2 [uniform, chord, or chord_square_root].
        Defaults to 0, uniform.

    Returns
    -------
    tuple (control_points, knot_vector)


    Raises
    ------
    ValueError
        If the knot style is not correct or if only one derivative is passed.
    """
    if knot_style not in [CurveKnotStyle.Uniform, CurveKnotStyle.Chord, CurveKnotStyle.ChordSquareRoot]:
        raise ValueError("Please pass a valid knot style: [0, 1, 2].")
    if start_derivative is not None and end_derivative is not None:
        return global_curve_interpolation_with_end_derivatives(points,
                                                               degree,
                                                               start_derivative,
                                                               end_derivative,
                                                               knot_style,
                                                               periodic=periodic)
    elif start_derivative is not None or end_derivative is not None:
        raise ValueError("Please pass start- AND end derivatives")
    else:
        return global_curve_interpolation(points, degree, knot_style=knot_style, periodic=periodic)


if __name__ == "__main__":
    from compas_nurbs import Curve

    knot_style = CurveKnotStyle.Uniform
    points = [[0.0, 0.0, 0.0], [0.973412, 1.962979, 0.0], [-0.66136, 2.434672, 0.0], [-2.34574, 1.016886, 0.0], [-3.513204, -0.855328, 0.0], [-4.0, -3.0, 0.0]]
    degree = 3

    control_points, knot_vector = global_curve_interpolation(points, degree, CurveKnotStyle.Uniform)
    crv = Curve(control_points, degree, knot_vector)

    start_derivative, end_derivative = [17.75, 10.79, 0.0], [-0.71, -12.62, 0.0]
    result = global_curve_interpolation_with_end_derivatives(points,
                                                             degree,
                                                             start_derivative,
                                                             end_derivative,
                                                             knot_style)

    crv = Curve(result[0], degree, result[1])
    D0, D1 = crv.derivatives_at([0, 1], order=1)[0]
    assert(allclose(D0, start_derivative))
    assert(allclose(D1, end_derivative))

    result = interpolate_curve(points, degree, knot_style)
    crv = Curve(result[0], degree, result[1])
    D0, D1 = crv.derivatives_at([0, 1], order=1)[0]
    assert(allclose(D0, start_derivative, tol=0.01))
    assert(allclose(D1, end_derivative, tol=0.01))
