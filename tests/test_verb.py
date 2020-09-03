from compas_nurbs import RationalSurface, Surface
from compas_nurbs import RationalCurve, Curve
from compas_nurbs.backends.verb import verb_core_NurbsSurfaceData
from compas_nurbs.backends.verb import verb_geom_NurbsSurface
from compas_nurbs.backends.verb import verb_geom_NurbsCurve
from compas_nurbs.backends.verb import verb_eval_Eval
from compas_nurbs.utilities import linspace
from compas.utilities import flatten
from compas.geometry import allclose, close


def test_nurbs_surface():
    # settings
    control_points_2d = [[[0, 0, 0, ], [0, 4, 0.], [0, 8, -3]],
                         [[2, 0, 6], [2, 4, 0.], [2, 8, 0.]],
                         [[4, 0, 0], [4, 4, 0.], [4, 8, 3.]],
                         [[6, 0, 0], [6, 4, -3], [6, 8, 0.]]]

    control_points_2d_homo = [[[x, y, z, 1.] for (x, y, z) in cl] for cl in control_points_2d]

    degree = (3, 2)
    weights = [[0.2, 0.1, 0.3], [0.1, 0.7, 1.2], [1.2, 2., 0.4], [0.1, 1.1, 0.5]]
    surface = RationalSurface(control_points_2d, degree, weights=weights)
    #surface = Surface(control_points_2d, degree)

    params_u = linspace(0., 1., 5)
    params_v = linspace(0., 1., 5)
    params = [(u, v) for u in params_u for v in params_v]

    degreeU, degreeV = surface.degree
    knotsU, knotsV = surface.knot_vector

    controlPoints = surface.control_points

    data = verb_core_NurbsSurfaceData(degreeU, degreeV, knotsU, knotsV, surface.weighted_control_points)
    verf_surface = verb_geom_NurbsSurface(data)

    p1 = [verb_eval_Eval.rationalSurfacePoint(data, u, v) for u, v in params]
    p2 = surface.points_at(params)
    assert(allclose(p1, p2))


def test_interpolation():
    pts = [[-10, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0], [5, 5, 0]]
    interpCurve = verb_geom_NurbsCurve.byPoints(pts, 3)
    crv = Curve.from_points(pts, 3, knot_style=1)  # chord style
    params = linspace(0., 1., 5)
    pts1 = [interpCurve.point(u) for u in params]
    pts2 = crv.points_at(params)
    assert(allclose(pts1, pts2))


def closest_point_to_curve():
    pts = [[-10, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0], [5, 5, 0]]
    interpCurve = verb_geom_NurbsCurve.byPoints(pts, 3)
    p1 = interpCurve.closestPoint([3., 5., 0.])
    t = interpCurve.closestParam(p1)
    p2 = interpCurve.point(t)
    assert(allclose(p1, p2))


if __name__ == "__main__":
    test_interpolation()
    closest_point_to_curve()
