import geomdl.BSpline
import geomdl.NURBS
from compas.utilities import flatten

# ==============================================================================
# curve
# ==============================================================================


def create_curve(control_points, degree, knot_vector, rational=False, weights=None):
    if not rational:
        curve = geomdl.BSpline.Curve()
        curve.degree = degree
        curve.ctrlpts = [list(pt) for pt in control_points]
        curve.knotvector = knot_vector
        return curve
    else:
        curve = geomdl.NURBS.Curve()
        curve.degree = degree
        curve.ctrlpts = [list(pt) for pt in control_points]
        curve.knotvector = knot_vector
        curve.weights = weights
        return curve


def evaluate_curve(curve, params):
    return curve.evaluate_list(params)


def evaluate_curve_derivatives(curve, params, order=1):
    derivatives = [curve.derivatives(u, order=order) for u in params]
    return derivatives

# ==============================================================================
# surface
# ==============================================================================


def create_surface(control_points_2d, degree, knot_vector, rational=False, weights=None):

    surface = geomdl.NURBS.Surface() if rational else geomdl.BSpline.Surface()
    surface.degree_u = degree[0]
    surface.degree_v = degree[1]
    count_u, count_v = len(control_points_2d), len(control_points_2d[0])
    ctrlpts = list(flatten(control_points_2d))

    if not rational:
        surface.set_ctrlpts(ctrlpts, count_u, count_v)
    else:
        weights = list(flatten(weights)) or [1. for _ in range(count_u, count_v)]
        ctrlptsw = [[w * x, w * y, w * z, w] for (x, y, z), w in zip(ctrlpts, weights)]
        surface.set_ctrlpts(ctrlptsw, count_u, count_v)

    surface.knotvector_u = knot_vector[0]
    surface.knotvector_v = knot_vector[1]

    return surface


def evaluate_surface(surface, params):
    return surface.evaluate_list(params)


def evaluate_surface_derivatives(surface, params, order=1):
    derivatives = [surface.derivatives(u, v, order=order) for (u, v) in params]
    return derivatives
