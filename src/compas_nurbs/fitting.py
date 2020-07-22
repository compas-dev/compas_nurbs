
# interpolation does not allow even degrees
# different rule with 3 [0000]

# 0.099897444571953278 * 10 = 1
# 17 / 5 = 0.29411764705882343 # degree 5 pts 13, cp 18
# 19 / 7 = 0.36842105263157893 # degree 7 pts 13, cp 20
# 21 / 9 = 0.42857142857142849 # degree 9 pts 13, cp 23


def compute_params_curve_equally_spaced(points):
    """ Computes :math:`\\overline{u}_{k}` for curves.

    Please refer to the equation 9.3 on The NURBS Book (2nd Edition), pp.364

    :param points: data points
    :type points: list, tuple
    :param centripetal: activates centripetal parametrization method
    :type centripetal: bool
    :return: parameter array, :math:`\\overline{u}_{k}`
    :rtype: list
    """
    n = len(points) - 1
    return [i / float(n) for i in range(len(points))]


def compute_knot_vector_with_end_derivatives(degree, num_points, params):
    """
    """
    kv = [0.0 for _ in range(degree + 1)]
    n = num_points - 1
    p = degree

    for j in range(n - p + 1 + 1):
        temp_kv = (1.0 / degree) * sum([params[j] for j in range(j, j + p - 1 + 1)])
        kv.append(temp_kv)

    kv += [1.0 for _ in range(degree + 1)]
    return kv


def global_curve_interpolation_with_end_tangents_numpy(points, degree,
                                                       knotstyle='uniform',
                                                       periodic=False,
                                                       start_tangent=None,
                                                       end_tangent=None):
    from scipy.linalg import lu_factor, lu_solve
    import numpy as np
    # TODO: does it anyway work with other tangents?
    start_tangent = start_tangent or [0, 0, 0]
    end_tangent = end_tangent or [0, 0, 0]

    num_points = len(points)
    uk = compute_params_curve_equally_spaced(points)
    kvd = compute_knot_vector_with_end_derivatives(degree, num_points, uk)  # isn't len(uk) == num_points
    M = _build_coeff_matrix_end_derivatives(degree, kvd, uk, points)
    C = points[:]
    C.insert(1, list(start_tangent))
    C.insert(-1, list(end_tangent))
    A = np.array(M)
    b = np.array(C)
    lu, piv = lu_factor(A)
    return lu_solve((lu, piv), b), kvd


def _build_coeff_matrix_end_derivatives(degree, knotvector, uk, points):
    """Builds the coefficient matrix for global interpolation with end derivatives specified.

    Please refer to The NURBS Book (2nd Edition), pp. 370-371 for details.

    Parameters
    ----------
    degree, int
    knotvector: list, tuple
    uk: list of parameters
    points: list, tuple

    Returns
    -------
    coefficient matrix, list
    """
    from geomdl import helpers
    num_points = len(points) + 2  # increase by 2

    # Set up coefficient matrix
    M = [[0.0 for _ in range(num_points)] for _ in range(num_points)]

    uke = [uk[0], 0.] + uk[1:-1] + [0., uk[-1]]  # fill with zeros

    for i in range(num_points):
        span = helpers.find_span_linear(degree, knotvector, num_points, uke[i])
        M[i][span - degree:span + 1] = helpers.basis_function(degree, knotvector, span, uke[i])

    M[1] = [0.0 for _ in range(num_points)]
    M[-2] = [0.0 for _ in range(num_points)]

    M[1][0] = -1.
    M[1][1] = 1.

    M[-2][-2] = -1.
    M[-2][-1] = 1.

    return M


if __name__ == "__main__":
    pass
