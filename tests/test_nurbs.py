from compas.geometry._core.average import weighted_centroid_points
import numpy as np
from compas.geometry import allclose


from compas_nurbs import NurbsCurve
from compas_nurbs import Curve
from compas_nurbs import knot_vectors
from compas_nurbs.utilities import linspace

from geomdl.NURBS import Curve as GNCurve
from geomdl.NURBS import Curve as GCurve


def test_bspline_curve():

    degree = 5
    control_points = [(0, 0, 0), (1.48, 0.9, 0), (4.95, 5.08, 0), (-1.56, 4.87, 0), (-4.72, -0.56, 0), (-3.94, -1.95, 0), (-4, -3, 0)]
    params = linspace(0., 1., 5)

    crv = Curve(control_points, degree)
    res1 = crv.points_at(params)

    ncrv = GCurve()
    ncrv.degree = degree
    ncrv.ctrlpts = [list(pt) for pt in control_points]
    ncrv.knotvector = crv.knot_vector

    res2 = ncrv.evaluate_list(params)
    assert(allclose(res1, res2))

    from scipy.interpolate import BSpline
    res3 = BSpline(crv.knot_vector, control_points, degree)(params)
    assert(allclose(res2, res3))


def test_nurbs_curve():
    degree = 5
    control_points = [(0, 0, 0), (1.48, 0.9, 0), (4.95, 5.08, 0), (-1.56, 4.87, 0), (-4.72, -0.56, 0), (-3.94, -1.95, 0), (-4, -3, 0)]
    weights = [0.5, 1.1, 0.7, 2., 4., 0.3, 0.1]
    params = linspace(0., 1., 5)

    crv = NurbsCurve(control_points, degree, weights=weights)
    res1 = crv.points_at(params)
    knot_vector = crv.knot_vector

    ncrv = GNCurve()
    ncrv.degree = degree
    ncrv.ctrlpts = [list(pt) for pt in control_points]
    ncrv.knotvector = knot_vector
    ncrv.weights = weights
    res2 = ncrv.evaluate_list(params)
    assert(allclose(res1, res2))

    d = [ncrv.derivatives(u, order=2)[1:] for u in params]
    print(d)
    D1 = list(zip(*d))
    D2 = crv.derivatives_at(params, order=2)
    assert(allclose(D1, D2))


if __name__ == "__main__":
    test_bspline_curve()
    test_nurbs_curve()
