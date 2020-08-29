import scipy
import numpy as np

from .helpers import find_spans
from .helpers import basis_functions

from compas.geometry._primitives.curve import binomial_coefficient  # TODO compas: move this upwards


def create_curve(control_points, degree, knot_vector, rational, weights):
    if not rational:
        return scipy.interpolate.BSpline(knot_vector, control_points, degree)
    else:
        w = np.array([weights]).T
        weighted_control_points = np.concatenate((w * control_points, w), axis=1)
        return scipy.interpolate.BSpline(knot_vector, weighted_control_points, degree)


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
        # TODO: make this faster with numpy
        derivatives.insert(0, evaluate_curve(curve, params, False))
        D = []
        for ders in zip(*derivatives):
            new_ders = []
            for k in range(0, order + 1):
                value = ders[k][:-1]
                for i in range(1, k + 1):
                    value[:] = [tmp - (binomial_coefficient(k, i) * ders[i][-1] * drv) for tmp, drv in zip(value, new_ders[k - i])]
                new_ders.append([tmp / ders[0][-1] for tmp in value])
            D.append(new_ders[1:])
        return np.array(D).transpose(1, 0, 2)


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
