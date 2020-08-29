import geomdl.BSpline
from compas.utilities import flatten


def create_curve(control_points, degree, knot_vector, rational, weights):
    if not rational:
        curve = geomdl.BSpline.Curve()
        curve.degree = degree
        curve.ctrlpts = [list(pt) for pt in control_points]
        curve.knotvector = knot_vector
        return curve
    else:
        curve = geomdl.Nurbs.Curve()
        curve.degree = degree
        curve.ctrlpts = [list(pt) for pt in control_points]
        curve.knotvector = knot_vector
        curve.weights = weights
        return curve


def evaluate_curve(curve, params, rational=False):
    return curve.evaluate_list(params)


def evaluate_curve_derivatives(curve, params, order=1):
    return curve.derivatives(params, order)
    # TODO


def create_surface(control_points_2d, degree_u, degree_v, knot_vector_u, knot_vector_v, rational=False, weights_u=None, weights_v=None):
    if not rational:
        surface = geomdl.BSpline.Surface()
        surface.degree_u = degree_u
        surface.degree_v = degree_v
        count_u, count_v = len(control_points_2d), len(control_points_2d[0])
        surface.set_ctrlpts(list(flatten(control_points_2d)), count_u, count_v)
        surface.knotvector_u = knot_vector_u
        surface.knotvector_v = knot_vector_v
        return surface
    else:
        pass


def evaluate_surface(surface, params, rational=False):
    return surface.evaluate_list(params)


def evaluate_surface_derivatives(surface, params, order=1):
    derivatives = [surface.derivatives(u, v, order=order) for (u, v) in params]
    D = []
    for i in range(1, order + 1):
        D.append([])
        for d in derivatives:
            D[i - 1].append(d[i])
    return D
