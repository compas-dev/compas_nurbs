import numpy as np
import rhino3dm
from geomdl import BSpline
from geomdl import NURBS
from compas.geometry import allclose
from compas.geometry import Point
from compas.geometry import Vector
from compas_nurbs import Curve
from compas_nurbs import RationalCurve


def rhino_curve_from_curve(curve):
    rhino_curve = rhino3dm.NurbsCurve(3, curve.rational, curve.degree + 1, curve.count)
    for i, p in enumerate(curve.control_points):
        rhino_curve.Points[i] = rhino3dm.Point4d(*p, curve.weights[i])
    for i, v in enumerate(curve.knot_vector[1:-1]):
        rhino_curve.Knots[i] = v
    return rhino_curve


def geomdl_curve_from_curve(curve):
    crv = NURBS.Curve() if curve.rational else BSpline.Curve()
    crv.degree = curve.degree
    crv.ctrlpts = curve.control_points
    crv.knotvector = curve.knot_vector
    crv.weights = curve.weights
    return crv


def test_curve():
    """Testing Curve against other libraries (geomdl, rhino3dm).
    """
    # settings
    control_points = [[816.32, 139.74, 0.0], [860.41, 238.68, 0.0], [928.17, 199.96, 0.0], [1008.83, 136.51, 0.0],
                      [1062.6, 195.66, 0.0], [1130.35, 248.36, 0.0], [1176.6, 202.11, 0.0], [1240.39, 164.93, 0.0],
                      [1315.67, 242.37, 0.0], [1396.33, 307.97, 0.0], [1466.23, 223.01, 0.0], [1520.01, 156.33, 0.0]]
    degree = 9
    params = np.linspace(0., 1., 10)

    # create curve
    curve = Curve(control_points, degree)
    curve_geomdl = geomdl_curve_from_curve(curve)
    curve_rhino = rhino_curve_from_curve(curve)

    # points
    points = curve.points_at(params)
    rhino_points = [curve_rhino.PointAt(t) for t in params]
    rhino_points = [[p.X, p.Y, p.Z] for p in rhino_points]
    geomdl_points = curve_geomdl.evaluate_list(params)
    assert(allclose(points, rhino_points))
    assert(allclose(points, geomdl_points))

    # derivatives
    geomdl_derivatives = [curve_geomdl.derivatives(u, order=1) for u in params]
    assert(allclose(geomdl_derivatives, curve.derivatives_at(params)))

    # tangents
    tangents = curve.tangents_at(params)
    rhino_tangents = [curve_rhino.TangentAt(t) for t in params]
    rhino_tangents = [[p.X, p.Y, p.Z] for p in rhino_tangents]
    geomdl_tangents = curve_geomdl.tangent(params)
    assert(allclose(tangents, rhino_tangents))
    assert(allclose(tangents, geomdl_tangents))

    # curvature
    curvature = [c.curvature for c in curve.curvatures_at(params)]
    rhino_curvature_vectors = [curve_rhino.CurvatureAt(t) for t in params]
    rhino_curvature = [np.linalg.norm([p.X, p.Y, p.Z]) for p in rhino_curvature_vectors]
    assert(allclose(curvature, rhino_curvature))
    circle_centers = [c.osculating_circle.plane.point for c in curve.curvatures_at(params)]
    rhino_centers = []
    for pt, cv in zip(rhino_points, rhino_curvature_vectors):
        cv = Vector(cv.X, cv.Y, cv.Z)
        rhino_centers.append(Point(*pt) + cv.unitized() * 1/cv.length)
    assert(allclose(circle_centers, rhino_centers))


def test_rational_curve():
    """Testing NurbsCurve against other libraries (geomdl, rhino3dm).

    Unfortunately the rhino3dm Nurbs Curve does not return the correct result.
    Instead rhino_* results for comparison are copied from Rhino directly.
    """
    # settings
    control_points = [[816.32, 139.74, 0.0], [860.41, 238.68, 0.0], [928.17, 199.96, 0.0], [1008.83, 136.51, 0.0],
                      [1062.6, 195.66, 0.0], [1130.35, 248.36, 0.0], [1176.6, 202.11, 0.0], [1240.39, 164.93, 0.0],
                      [1315.67, 242.37, 0.0], [1396.33, 307.97, 0.0], [1466.23, 223.01, 0.0], [1520.01, 156.33, 0.0]]
    degree = 9
    params = np.linspace(0., 1., 10)
    weights = [0.5, 1.1, 0.7, 2., 4., 0.3, 0.1, 0.3, 0.4, 0.5, 1.4, 3.]

    # create curve
    curve = RationalCurve(control_points, degree, weights=weights)
    curve_geomdl = geomdl_curve_from_curve(curve)
    # curve_rhino = rhino_curve_from_curve(curve)

    # points
    points = curve.points_at(params)
    rhino_points = [(816.320, 139.740, 0.), (944.679, 193.328, 0.), (1020.014, 175.383, 0.), (1045.948, 181.534, 0.),
                    (1065.127, 189.409, 0.), (1109.071, 198.051, 0.), (1226.082, 218.523, 0.), (1348.823, 250.108, 0.),
                    (1438.406, 239.541, 0.), (1520.010, 156.330, 0.)]
    geomdl_points = curve_geomdl.evaluate_list(params)
    assert(allclose(points, rhino_points, tol=1e-03))
    assert(allclose(points, geomdl_points))

    # derivatives
    geomdl_derivatives = [curve_geomdl.derivatives(u, order=1) for u in params]
    assert(allclose(geomdl_derivatives, curve.derivatives_at(params)))

    # tangents
    tangents = curve.tangents_at(params)
    rhino_tangents = [(0.407, 0.913, 0.), (0.939, -0.345, 0.), (0.999, 0.033, 0.), (0.917, 0.400, 0.),
                      (0.954, 0.301, 0.), (0.989, 0.151, 0.), (0.977, 0.213, 0.), (0.977, 0.212, 0.),
                      (0.855, -0.518, 0.), (0.628, -0.778, 0.)]
    geomdl_tangents = curve_geomdl.tangent(params)
    assert(allclose(tangents, rhino_tangents, tol=1e-03))
    assert(allclose(tangents, geomdl_tangents))

    # curvature
    curvature = [c.curvature for c in curve.curvatures_at(params)]
    rhino_curvature = [0.000851, 0.000762, 0.013305, 0.006902, 0.009795, 0.00037, 0.000876, 0.003371, 0.008178, 5e-05]
    assert(allclose(curvature, rhino_curvature))
    circle_centers = [c.osculating_circle.plane.point for c in curve.curvatures_at(params)]
    rhino_centers = [(1889.858, -338.654, 0.), (1396.986, 1424.413, 0.), (1017.561, 250.506, 0.), (988.023, 314.332, 0.),
                     (1095.819, 92.043, 0.), (1516.222, -2474.241, 0.), (982.641, 1333.837, 0.), (1411.820, -39.776, 0.),
                     (1375.058, 134.945, 0.), (-14154.272, -12485.585, 0.)]
    assert(allclose(circle_centers, rhino_centers, tol=1e-03))


if __name__ == "__main__":
    test_curve()
    test_rational_curve()
