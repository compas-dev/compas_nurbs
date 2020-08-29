import numpy as np
import rhino3dm
from compas.geometry import allclose
from compas_nurbs import Curve
from compas_nurbs import Surface
from compas_nurbs.knot_vectors import knot_vector_uniform


def test_compare_bspline_with_rhino_curve():
    """The default rhino curve evaluates the same way as BSpline.from_uniform_knot_style.
    """

    control_points = [[816.32, 139.74, 0.0], [860.41, 238.68, 0.0], [928.17, 199.96, 0.0], [1008.83, 136.51, 0.0],
                      [1062.6, 195.66, 0.0], [1130.35, 248.36, 0.0], [1176.6, 202.11, 0.0], [1240.39, 164.93, 0.0],
                      [1315.67, 242.37, 0.0], [1396.33, 307.97, 0.0], [1466.23, 223.01, 0.0], [1520.01, 156.33, 0.0]]
    degree = 9
    #control_points = [(0.68, 0.44, 0.00), (0.27, 2.50, 0.00), (6.03, 2.18, 0.00), (4.77, 4.50, 0.00)]
    #degree = 3

    n = 10
    params = [i / float(n) for i in range(n + 1)]

    # create BSpline
    curve = Curve(control_points, degree)

    # create Rhino curve
    rhino_control_points = [rhino3dm.Point3d(*p) for p in control_points]
    rhino_curve = rhino3dm.Curve.CreateControlPointCurve(rhino_control_points, degree)
    rhino_params = [t * rhino_curve.Domain.T1 for t in params]  # rhino curve domain is different

    # evaluate points
    points = curve.points_at(params)
    rhino_points = [rhino_curve.PointAt(t) for t in rhino_params]
    rhino_points = [[p.X, p.Y, p.Z] for p in rhino_points]

    # evaluate tangents
    tangents = curve.tangents_at(params)
    rhino_tangents = [rhino_curve.TangentAt(t) for t in rhino_params]
    rhino_tangents = [[p.X, p.Y, p.Z] for p in rhino_tangents]

    # evaluate curvature
    curvature = curve.curvatures_at(params)
    rhino_curvature = [rhino_curve.CurvatureAt(t) for t in rhino_params]
    rhino_curvature = [np.linalg.norm([p.X, p.Y, p.Z]) for p in rhino_curvature]

    # compare
    assert(allclose(points, rhino_points))
    assert(allclose(tangents, rhino_tangents))
    assert(allclose(curvature, rhino_curvature))

    # compare with scipy
    from scipy.interpolate import BSpline
    scipy_curve = BSpline(curve.knot_vector, curve.control_points, curve.degree)
    scipy_points = scipy_curve(params)

    print(tangents)
    d1 = scipy_curve.derivative(1)(params)
    scipy_tangents = (d1.T/np.linalg.norm(d1, axis=1).reshape(1, -1)).T

    d2 = scipy_curve.derivative(2)(params)
    k = np.linalg.norm(np.cross(d1, d2, axis=1), axis=1) / np.linalg.norm(d1, axis=1)**3
    scipy_curvature = k

    assert(allclose(points, scipy_points))
    assert(allclose(tangents, scipy_tangents))
    assert(allclose(curvature, scipy_curvature))


"""
def test_compare_surface_with_geomdl_surface():
    from compas.geometry import allclose
    from compas_nurbs.utilities import linspace

    degree_u, degree_v = 3, 3
    count_u, count_v = 6, 6
    control_points = [[-25.0, -25.0, -10.0], [-25.0, -15.0, -5.0], [-25.0, -5.0, 0.0], [-25.0, 5.0, 0.0], [-25.0, 15.0, -5.0], [-25.0, 25.0, -10.0],
                      [-15.0, -25.0, -8.0], [-15.0, -15.0, -4.0], [-15.0, -5.0, -4.0], [-15.0, 5.0, -4.0], [-15.0, 15.0, -4.0], [-15.0, 25.0, -8.0],
                      [-5.0, -25.0, -5.0], [-5.0, -15.0, -3.0], [-5.0, -5.0, -8.0], [-5.0, 5.0, -8.0], [-5.0, 15.0, -3.0], [-5.0, 25.0, -5.0],
                      [5.0, -25.0, -3.0], [5.0, -15.0, -2.0], [5.0, -5.0, -8.0], [5.0, 5.0, -8.0], [5.0, 15.0, -2.0], [5.0, 25.0, -3.0],
                      [15.0, -25.0, -8.0], [15.0, -15.0, -4.0], [15.0, -5.0, -4.0], [15.0, 5.0, -4.0], [15.0, 15.0, -4.0], [15.0, 25.0, -8.0],
                      [25.0, -25.0, -10.0], [25.0, -15.0, -5.0], [25.0, -5.0, 2.0], [25.0, 5.0, 2.0], [25.0, 15.0, -5.0], [25.0, 25.0, -10.0]]
    knot_vector_u = knot_vector_uniform(degree_u, count_u)
    knot_vector_v = knot_vector_uniform(degree_v, count_v)
    knot_vector_v = knot_vector_u
    srf = Surface(control_points, count_u, count_v, degree_u, degree_v, knot_vector_u, knot_vector_v)

    params_u = linspace(0., 1., 5)
    params_v = linspace(0., 1., 5)
    params_uv = [(u, v) for u in params_u for v in params_v]

    points = srf.evaluate_at(params_u, params_v)

    from geomdl import BSpline
    surf = BSpline.Surface()
    surf.degree_u = degree_u
    surf.degree_v = degree_v
    surf.set_ctrlpts(control_points, count_u, count_v)
    surf.knotvector_u = knot_vector_u
    surf.knotvector_v = knot_vector_v

    geomdl_points = surf.evaluate_list(params_uv)
    assert(allclose(points, geomdl_points))
"""


def test_compare_interpolation_with_scipy_interpolation():
    """Scipy bspline interpolation only till degree 5
    """
    import numpy as np
    import scipy.interpolate as interpolate
    import matplotlib.pyplot as plt

    x = np.array([0.,  1.2,  1.9,  3.2,  4.,  6.5])
    y = np.array([0.,  2.3,  3.,  4.3,  2.9,  3.1])

    t, c, k = interpolate.splrep(x, y, s=0, k=4)
    # t = vector of knots
    # c = the B-spline coefficients
    # k = degree of the spline

    N = 100
    xmin, xmax = x.min(), x.max()
    xx = np.linspace(xmin, xmax, N)
    spline = interpolate.BSpline(t, c, k, extrapolate=False)
    # plt.plot(xx, spline(xx)
    # TODO


def test_compare_evaluators():
    import os
    import json
    from compas_nurbs import DATA
    from compas_nurbs.evaluators import evaluate_curve
    from compas_nurbs.evaluators_numpy import evaluate_curve as evaluate_curve_numpy
    import geomdl.BSpline

    with open(os.path.join(DATA, "bsplines_data.json"), 'r') as f:
        data = json.load(f)

    curve = Curve.from_data(data)
    params = np.linspace(0, 1, 1000)

    geomdl_curve = geomdl.BSpline.Curve()
    geomdl_curve.degree = curve.degree
    geomdl_curve.ctrlpts = [list(pt) for pt in curve.control_points]
    geomdl_curve.knotvector = curve.knot_vector

    pts1 = evaluate_curve(geomdl_curve, params, rational=False)
    pts2 = evaluate_curve_numpy(curve._curve, params, rational=False)
    assert(allclose(pts1, pts2))


if __name__ == "__main__":
    test_compare_bspline_with_rhino_curve()
    test_compare_evaluators()
