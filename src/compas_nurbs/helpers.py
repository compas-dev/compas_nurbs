def find_span(knot_vector, number_of_control_points, t):
    """Finds the span of a single knot over the knot vector using linear search.

    Alternative implementation for the Algorithm A2.1 from The NURBS Book by Piegl & Tiller.

    # TODO find_span_binsearch
    """
    span = 0
    while span < number_of_control_points and knot_vector[span] <= t:
        span += 1
    return span - 1


def find_spans(knot_vector, number_of_control_points, params):
    """
    TODO: check if np.vectorize makes it faster.
    """
    return list(map(lambda t: find_span(knot_vector, number_of_control_points, t), params))


def basis_function(degree, knot_vector, span, t):
    """Computes the non-vanishing basis functions for a single parameter t.

    Implementation of Algorithm A2.2 from The NURBS Book by Piegl & Tiller.
    Uses recurrence to compute the basis functions, also known as Cox - de
    Boor recursion formula.
    """
    left = [0.0 for _ in range(degree + 1)]
    right = [0.0 for _ in range(degree + 1)]
    N = [1.0 for _ in range(degree + 1)]

    for j in range(1, degree + 1):
        left[j] = t - knot_vector[span + 1 - j]
        right[j] = knot_vector[span + j] - t
        saved = 0.0
        for r in range(0, j):
            temp = N[r] / (right[r + 1] + left[j - r])
            N[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        N[j] = saved
    return N


def basis_function_derivatives(degree, knot_vector, span, t, order):
    """Computes derivatives of the basis functions for a single parameter t.

    Implementation of Algorithm A2.3 from The NURBS Book by Piegl & Tiller.

    Parameters
    ----------
    degree: int
    knot_vector: knot vector
    span: int
    t: float
    order: int
        The order of the derivative.

    Returns
    -------
    list : derivatives of the basis functions
    """
    # Initialize variables
    left = [1.0 for _ in range(degree + 1)]
    right = [1.0 for _ in range(degree + 1)]
    ndu = [[1.0 for _ in range(degree + 1)] for _ in range(degree + 1)]  # N[0][0] = 1.0 by definition

    for j in range(1, degree + 1):
        left[j] = t - knot_vector[span + 1 - j]
        right[j] = knot_vector[span + j] - t
        saved = 0.0
        r = 0
        for r in range(r, j):
            # Lower triangle
            ndu[j][r] = right[r + 1] + left[j - r]
            temp = ndu[r][j - 1] / ndu[j][r]
            # Upper triangle
            ndu[r][j] = saved + (right[r + 1] * temp)
            saved = left[j - r] * temp
        ndu[j][j] = saved

    # Load the basis functions
    ders = [[0.0 for _ in range(degree + 1)] for _ in range((min(degree, order) + 1))]
    for j in range(0, degree + 1):
        ders[0][j] = ndu[j][degree]

    # Start calculating derivatives
    a = [[1.0 for _ in range(degree + 1)] for _ in range(2)]
    # Loop over function index
    for r in range(0, degree + 1):
        # Alternate rows in array a
        s1 = 0
        s2 = 1
        a[0][0] = 1.0
        # Loop to compute k-th derivative
        for k in range(1, order + 1):
            d = 0.0
            rk = r - k
            pk = degree - k
            if r >= k:
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk]
                d = a[s2][0] * ndu[rk][pk]
            if rk >= -1:
                j1 = 1
            else:
                j1 = -rk
            if (r - 1) <= pk:
                j2 = k - 1
            else:
                j2 = degree - r
            for j in range(j1, j2 + 1):
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j]
                d += (a[s2][j] * ndu[rk + j][pk])
            if r <= pk:
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r]
                d += (a[s2][k] * ndu[r][pk])
            ders[k][r] = d

            # Switch rows
            j = s1
            s1 = s2
            s2 = j

    r = float(degree)  # Multiply through by the the correct factors
    for k in range(1, order + 1):
        for j in range(0, degree + 1):
            ders[k][j] *= r
        r *= (degree - k)
    return ders


def basis_functions(degree, knot_vector, spans, params):
    """
    TODO: check if np.vectorize makes it faster.
    """
    return list(map(lambda span, t: basis_function(degree, knot_vector, span, t), spans, params))


def basis_functions_derivatives(degree, knot_vector, spans, params, order):
    """
    TODO: check if np.vectorize makes it faster.
    """
    return list(map(lambda span, t: basis_function_derivatives(degree, knot_vector, span, t, order), spans, params))


"""

or base on scipy's bspline?
https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.BSpline.html

BSline.basis_element(t[, extrapolate])activa
BSline.derivative(self[, nu])
BSline.integrate(self, a, b[, extrapolate])

"""
