import numpy as np

from .helpers import find_spans
from .helpers import basis_functions
from .helpers import basis_functions_derivatives


def evaluate_curve(curve, params, rational=False):
    """Evaluates a curve at the parameters params.
    """
    points = curve(params)
    
    if not rational:
        return points
    else:
        w = np.array([points[:, -1]]).T
        return np.delete((1/w * points), -1, axis=1)


def evaluate_curve_derivatives(curve, params, order=1, rational=False):
    """Evaluates the n-th order derivatives at the parametric positions `params`.

    Parameters
    ----------
    control_points: np.array
    degree : int
        The degree of the curve
    knot_vector : list
        The knot vector of the curve
    params : list
        Parametric positions where the derivatives will be computed. Range [0, 1]
    order : int
        The derivative order; to get the i-th derivative

    Returns
    -------
    :class:`numpy.array`
        The evaluated derivatives

    Examples
    --------
    >>>
    """
    derivatives = []
    for i in range(1, order + 1):
        derivatives.append(curve.derivative(i)(params))
    
    if not rational:
        return derivatives
    else:
        print(derivatives)
        print("NotImplementedError")
        #w = np.array([derivatives[:, -1]]).T
        #return np.delete((1/w * derivatives), -1, axis=1)

    """
    # Call the parent function to evaluate A(u) and w(u) derivatives
        CKw = super(CurveEvaluatorRational, self).derivatives(datadict, parpos, deriv_order, **kwargs)

        # Algorithm A4.2
        CK = [[0.0 for _ in range(dimension - 1)] for _ in range(deriv_order + 1)]
        for k in range(0, deriv_order + 1):
            v = [val for val in CKw[k][0:(dimension - 1)]]
            for i in range(1, k + 1):
                v[:] = [tmp - (linalg.binomial_coefficient(k, i) * CKw[i][-1] * drv) for tmp, drv in
                        zip(v, CK[k - i])]
            CK[k][:] = [tmp / CKw[0][-1] for tmp in v]

        # Return C(u) derivatives
        return CK
    """


def evaluate_surface(control_points, degree_u, degree_v, knot_vector_u, knot_vector_v, params_u, params_v, rational=False):
    """Evaluates a surface at the parameters params_u, params_v.
    """
    if rational:
        raise NotImplementedError

    number_of_control_points_u = len(control_points)
    number_of_control_points_v = len(control_points[0])
    spans_u = find_spans(knot_vector_u, number_of_control_points_u, params_u)
    bases_u = basis_functions(degree_u, knot_vector_u, spans_u, params_u)
    spans_v = find_spans(knot_vector_v, number_of_control_points_v, params_v)
    bases_v = basis_functions(degree_v, knot_vector_v, spans_v, params_v)

    points = []
    for span_u, basis_u in zip(spans_u, bases_u):  # any chance to make this faster with matrix multiplication?
        idx_u = span_u - degree_u
        for span_v, basis_v in zip(spans_v, bases_v):
            idx_v = span_v - degree_v
            a = basis_v[:degree_v + 1]
            b = control_points[idx_u:idx_u + degree_u + 1, idx_v:idx_v + degree_v + 1]
            c = basis_u[:degree_u + 1]
            points.append(np.dot(c, np.dot(a, b)))

    return np.array(points)


def evaluate_surface_derivatives(control_points, degree_u, degree_v, knot_vector_u, knot_vector_v, params_u, params_v, order=1, normalize=True):
    raise NotImplementedError
