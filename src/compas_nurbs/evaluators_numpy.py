import numpy as np
from .helpers import find_spans
from .helpers import basis_functions
from .helpers import basis_functions_derivatives


def evaluate_curve(control_points, degree, knot_vector, params, rational=False):
    """Evaluates a curve at the parameters params.
    """
    if rational:
        raise NotImplementedError

    number_of_control_points = len(control_points)
    spans = find_spans(knot_vector, number_of_control_points, params)
    bases = np.array(basis_functions(degree, knot_vector, spans, params))
    control_points = np.array(control_points)

    points = []
    for span, basis in zip(spans, bases):  # any chance to make this faster with matrix multiplication?
        a = basis[:degree + 1]
        b = control_points[span - degree:span + 1]
        points.append(np.dot(a, b))
    return np.array(points)


def evaluate_curve_derivatives(control_points, degree, knot_vector, params, order=1, normalize=True):
    """Evaluates the n-th order derivatives at the input parametric position.

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
    number_of_control_points = len(control_points)

    if order > degree:
        raise ValueError("Order must be < than degree.")

    spans = find_spans(knot_vector, number_of_control_points, params)
    bases = basis_functions_derivatives(degree, knot_vector, spans, params, order)
    bases = np.array(bases)[:, order, :]

    derivatives = []
    for span, basis in zip(spans, bases):  # any chance to make this faster with matrix multiplication?
        a = basis[:degree + 1]
        b = control_points[span - degree:span + 1]
        derivatives.append(np.dot(a, b))
    derivatives = np.array(derivatives)
    if normalize:
        return (derivatives.T / np.linalg.norm(derivatives, axis=1)).T
    else:
        return derivatives


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
    pass
