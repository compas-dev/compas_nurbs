import scipy
import numpy as np

from .helpers import find_spans
from .helpers import basis_functions
from .helpers import basis_functions_derivatives

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
    for i in range(0, order + 1):
        derivatives.append(curve.derivative(i)(params))
    derivatives = np.array(derivatives).transpose(1, 0, 2)

    if not rational:
        return derivatives
    else:
        # TODO: numpify this!
        #derivatives.insert(0, evaluate_curve(curve, params, False))
        D = []
        # for ders in zip(*derivatives):
        for ders in derivatives:
            new_ders = []
            for k in range(0, order + 1):
                value = ders[k][:-1]
                for i in range(1, k + 1):
                    value[:] = [tmp - (binomial_coefficient(k, i) * ders[i][-1] * drv) for tmp, drv in zip(value, new_ders[k - i])]
                new_ders.append([tmp / ders[0][-1] for tmp in value])
            D.append(new_ders[1:])
        return np.array(D).transpose(1, 0, 2)

# ===============================================================================
# surface
# ===============================================================================


def create_surface(control_points_2d, degree_u, degree_v, knot_vector_u, knot_vector_v, rational, weights_u, weights_v):
    if not rational:
        pass
    else:
        pass


def evaluate_surface(surface, params, rational=False):
    """Evaluates a surface at the parameters params_u, params_v.
    """
    if rational:
        raise NotImplementedError

    control_points = np.array(surface.control_points_2d)
    degree_u, degree_v = surface.degree_u, surface.degree_v
    knot_vector_u, knot_vector_v = surface.knot_vector_u, surface.knot_vector_v

    number_of_control_points_u = len(control_points)
    number_of_control_points_v = len(control_points[0])

    params_u = [p[0] for p in params]
    params_v = [p[1] for p in params]

    spans_u = find_spans(knot_vector_u, number_of_control_points_u, params_u)
    bases_u = basis_functions(degree_u, knot_vector_u, spans_u, params_u)
    spans_v = find_spans(knot_vector_v, number_of_control_points_v, params_v)
    bases_v = basis_functions(degree_v, knot_vector_v, spans_v, params_v)

    points = []
    for span_u, basis_u, span_v, basis_v in zip(spans_u, bases_u, spans_v, bases_v):
        a = basis_v[:degree_v + 1]
        b = control_points[span_u - degree_u:span_u + 1, span_v - degree_v:span_v+1]
        c = basis_u[:degree_u + 1]
        points.append(np.dot(c, np.dot(a, b)))
    return np.array(points)


def evaluate_surface_derivatives(surface, params, order=1):
    """
    """
    control_points = np.array(surface.control_points_2d)
    degree_u, degree_v = surface.degree_u, surface.degree_v
    knot_vector_u, knot_vector_v = surface.knot_vector_u, surface.knot_vector_v

    number_of_control_points_u = len(control_points)
    number_of_control_points_v = len(control_points[0])

    params_u = [p[0] for p in params]
    params_v = [p[1] for p in params]

    spans_u = find_spans(knot_vector_u, number_of_control_points_u, params_u)
    bases_u = basis_functions_derivatives(degree_u, knot_vector_u, spans_u, params_u, order)
    spans_v = find_spans(knot_vector_v, number_of_control_points_v, params_v)
    bases_v = basis_functions_derivatives(degree_v, knot_vector_v, spans_v, params_v, order)

    dv = min(degree_v, order)
    derivatives = []
    for span_u, basis_u, span_v, basis_v in zip(spans_u, bases_u, spans_v, bases_v):
        b = control_points[span_u - degree_u:span_u + 1, span_v - degree_v:span_v + 1]
        temp = np.dot(b.T, np.array(basis_u).T).T
        dd = min(order, dv)
        SKL = np.dot(np.array(basis_v[:dd+1]), temp[:degree_v + 1]).transpose(1, 0, 2)
        derivatives.append(SKL)
    return np.array(derivatives)


def calculate_surface_curvature(derivatives, order=False):
    """Calculates surface curvature quantities.

    If "order" parameter is set to True, then it will be guaranteed, that C1
    value is always less than C2.
    """
    fu = derivatives[:, 1, 0]
    fv = derivatives[:, 0, 1]
    fuu = derivatives[:, 2, 0]
    fvv = derivatives[:, 0, 2]
    fuv = derivatives[:, 1, 1]

    normal = np.cross(fu, fv)
    norm = np.linalg.norm(normal, axis=1, keepdims=True)
    normal = normal / norm

    nuu = (fuu * normal).sum(axis=1)
    nvv = (fvv * normal).sum(axis=1)
    nuv = (fuv * normal).sum(axis=1)

    duu = np.linalg.norm(fu, axis=1) ** 2
    dvv = np.linalg.norm(fv, axis=1) ** 2
    duv = (fu * fv).sum(axis=1)

    mean = (duu*nvv - 2*duv*nuv + dvv*nuu) / (2*(duu*dvv - duv*duv))
    gauss = (nuu * nvv - nuv*nuv) / (duu * dvv - duv*duv)

    n = len(derivatives)  # number of params
    L = np.empty((n, 2, 2))
    L[:, 0, 0] = nuu
    L[:, 0, 1] = nuv
    L[:, 1, 0] = nuv
    L[:, 1, 1] = nvv

    G = np.empty((n, 2, 2))
    G[:, 0, 0] = duu
    G[:, 0, 1] = duv
    G[:, 1, 0] = duv
    G[:, 1, 1] = dvv

    M = np.matmul(np.linalg.inv(G), L)
    eigvals, eigvecs = np.linalg.eig(M)
    # Values of first and second principal curvatures
    c1 = eigvals[:, 0]
    c2 = eigvals[:, 1]

    if order:
        c1mask = (c1 < c2)
        c2mask = np.logical_not(c1mask)
        c1_r = np.where(c1mask, c1, c2)
        c2_r = np.where(c2mask, c1, c2)
    else:
        c1_r = c1
        c2_r = c2

    # dir_1 corresponds to c1, dir_2 corresponds to c2
    dir_1_x = eigvecs[:, 0, 0][np.newaxis].T
    dir_2_x = eigvecs[:, 0, 1][np.newaxis].T
    dir_1_y = eigvecs[:, 1, 0][np.newaxis].T
    dir_2_y = eigvecs[:, 1, 1][np.newaxis].T
    dir_1 = dir_1_x * fu + dir_1_y * fv
    dir_2 = dir_2_x * fu + dir_2_y * fv
    dir_1 = dir_1 / np.linalg.norm(dir_1, axis=1, keepdims=True)
    dir_2 = dir_2 / np.linalg.norm(dir_2, axis=1, keepdims=True)

    if order:
        c1maskT = c1mask[np.newaxis].T
        c2maskT = c2mask[np.newaxis].T
        dir_1_r = np.where(c1maskT, dir_1, -dir_2)
        dir_2_r = np.where(c2maskT, dir_1, dir_2)
    else:
        dir_1_r = dir_1
        dir_2_r = dir_2

    return c1_r, c2_r, dir_1_r, dir_2_r, normal, mean, gauss
