def evaluate_curve(curve, params, rational=False):
    return curve.evaluate_list(params)


def evaluate_curve_derivatives(curve, params, order=1):
    return curve.derivatives(params, order)


def evaluate_surface(control_points, degree_u, degree_v, knot_vector_u, knot_vector_v, params_u, params_v, rational=False):
    from geomdl.BSpline import Surface as GSurface
    srf = GSurface()
    srf.degree_u = degree_u
    srf.degree_v = degree_v
    nu = len(control_points)
    nv = len(control_points[0])
    control_points_flattened = [list(pt) for pts in control_points for pt in pts]
    srf.set_ctrlpts(control_points_flattened, nu, nv)
    srf.knotvector_u = knot_vector_u
    srf.knotvector_v = knot_vector_v
    params = [(u, v) for v in params_v for u in params_u]
    return srf.evaluate_list(params)


def evaluate_surface_derivatives(control_points, degree_u, degree_v, knot_vector_u, knot_vector_v, params_u, params_v, order=1, normalize=True):
    raise NotImplementedError
