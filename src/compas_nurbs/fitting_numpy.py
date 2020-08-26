import numpy as np
from scipy.linalg import lu_factor, lu_solve

from compas_nurbs.helpers import find_spans
from compas_nurbs.helpers import basis_functions
from compas_nurbs.knot_vectors import knot_vector_and_params



# https://link.springer.com/content/pdf/10.1007/s003660050038.pdf



def global_curve_interpolation(points, degree, knot_style='uniform', periodic=False):
    """Curve interpolation through points.

    Please refer to Algorithm A9.1 on The NURBS Book (2nd Edition), pp.369-370
    for details.
    
    Parameters
    ----------
    points : list of point
    degree : int
        degree of the output parametric curve
    knotstyle : str, optional
        Defaults to 'uniform'
    
    Returns
    -------
    tuple (control_points, knot_vector)
        The interpolated B-Spline curve
    

    Examples
    --------
    >>> degree = 3
    >>> points = [(0., 0., 0.), (3., 4., 0.), (-1., 4., 0.), (-4., 0., 0.), (-4., -3., 0.)]
    >>> control_points, knot_vector = global_curve_interpolation(points, degree)
    >>> solution = [(0, 0, 0), (7.32, 3.69, 0), (-2.96, 6.68, 0), (-4.49, -0.67, 0), (-4, -3, 0)]
    >>> allclose(control_points, solution, tol=0.01)
    True
    """
    points = np.array(points)

    kv, uk = knot_vector_and_params(points, degree, knot_style, extended=False)

    M = coefficient_matrix(degree, kv, uk)
    lu, piv = lu_factor(M)
    return lu_solve((lu, piv), points), kv


def estimate_end_derivates():
    # create a curve where 1st and 2nd point are the same
    pass


"""
If unit length tangent vectors Tk are given instead of the Dk, then magnitudes
uk must be estimated. Setting all ak = d, where d is the total chord length, is
a reasonable choice.
"""






def global_curve_interpolation_with_end_derivatives(points, degree,
                                                       knot_style='uniform',
                                                       periodic=False,
                                                       start_derivative=None,
                                                       end_derivative=None):

    points = np.array(points)

    kv, uk = knot_vector_and_params(points, degree, knot_style, extended=True)

    # Please refer to The NURBS Book (2nd Edition), pp. 370-371 for details.
    M = coefficient_matrix(degree, kv, uk)
    M[1][0] = -1.
    M[1][1] = 1.
    M[-2][-2] = -1.
    M[-2][-1] = 1.

    D10 = start_derivative
    D1n = end_derivative
    v0 = np.array(D10) * kv[degree+1]/degree
    vn = np.array(D1n) * (1 - kv[len(kv) - 1 - degree-1])/degree
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
    knot_vector : list of float
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


def interpolate_curve(points, degree, knot_style='uniform', start_derivative=None, end_derivative=None, periodic=False):
    """
    """

    if degree == 3:
        # estimate end tangents
        pass
    elif start_derivative and end_derivative:
        return global_curve_interpolation_with_end_derivatives(points, degree,
                                                       knot_style=knot_style,
                                                       periodic=periodic,
                                                       start_derivative=start_derivative,
                                                       end_derivative=end_derivative)
    else:
        return global_curve_interpolation(points, degree, knot_style=knot_style, periodic=periodic)



if __name__ == "__main__":
    from compas.geometry import allclose # noqa: F401
    points = [[0.0, 0.0, 0.0], [0.973412, 1.962979, 0.0], [-0.66136, 2.434672, 0.0], [-2.34574, 1.016886, 0.0], [-3.513204, -0.855328, 0.0], [-4.0, -3.0, 0.0]]
    degree = 7
    cpts, kv = global_curve_interpolation_with_end_tangents_numpy(points, degree, )
    solution = [[0.0, 0.0, 0.0], [3.855235, 0.75935299999999994, 0.0], [-0.053739000000000002, 4.2156969999999996, 0.0], [-2.4180000000000001, 3.8224860000000001, 0.0],
              [-2.3438430000000001, 0.73678500000000002, 0.0], [-2.5536300000000001, -1.20214, 0.0], [-4.5866569999999998, -0.49421999999999999, 0.0], [-4.0, -3.0, 0.0]]
    print(allclose(cpts, solution))

    from compas_nurbs import Curve
    
    curve = Curve.from_uniform_knot_style(solution, degree=degree)
    print(curve.derivatives_at([0, 1], order=1))

    """
    [[ 0.          0.          0.        ]
    [ 3.855235    0.759353    0.        ]
    [-1.69592187  4.77690428  0.        ]
    [ 1.52264906  2.74022362  0.        ]
    [-4.89288482  2.14450472  0.        ]
    [-1.69889489 -1.41738231  0.        ]
    [-4.586657   -0.49422    -0.        ]
    [-4.         -3.          0.        ]]
    """

    print(curve.knot_vector)
    


    #import doctest
    #from compas.geometry import allclose  
    #doctest.testmod(globs=globals())