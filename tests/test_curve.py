import numpy as np
import rhino3dm
from compas.geometry import allclose
from compas_nurbs import Curve
from compas_nurbs import NurbsCurve
from compas_nurbs.evaluators import create_curve
from compas_nurbs.evaluators import evaluate_curve
from compas_nurbs.evaluators import evaluate_curve_derivatives
from compas_nurbs.operations import curve_tangents
from compas_nurbs.utilities import linspace


def rhino_curve_from_curve(curve, rational=False):
    rhino_curve = rhino3dm.NurbsCurve(3, rational, curve.degree + 1, curve.count)
    for i, p in enumerate(curve.control_points):
        x, y, z = p
        if rational:
            rhino_curve.Points[i] = rhino3dm.Point4d(x, y, z, curve.weights[i])
        else:
            rhino_curve.Points[i] = rhino3dm.Point4d(x, y, z, 1.)

    for i, v in enumerate(curve.knot_vector[1:-1]):
        rhino_curve.Knots[i] = v
    return rhino_curve


def test_curve():
    """Testing Curve against other libraries (geomdl, rhino3dm).
    """
    # settings
    control_points = [[816.32, 139.74, 0.0], [860.41, 238.68, 0.0], [928.17, 199.96, 0.0], [1008.83, 136.51, 0.0],
                      [1062.6, 195.66, 0.0], [1130.35, 248.36, 0.0], [1176.6, 202.11, 0.0], [1240.39, 164.93, 0.0],
                      [1315.67, 242.37, 0.0], [1396.33, 307.97, 0.0], [1466.23, 223.01, 0.0], [1520.01, 156.33, 0.0]]
    degree = 9
    params = linspace(0., 1., 10)

    # create curve
    curve = Curve(control_points, degree)
    curve_geomdl = create_curve(control_points, degree, curve.knot_vector)
    curve_rhino = rhino_curve_from_curve(curve)

    # points
    points = curve.points_at(params)
    rhino_points = [curve_rhino.PointAt(t) for t in params]
    rhino_points = [[p.X, p.Y, p.Z] for p in rhino_points]
    geomdl_points = evaluate_curve(curve_geomdl, params)
    assert(allclose(points, rhino_points))
    assert(allclose(points, geomdl_points))

    # derivatives
    geomdl_derivatives = evaluate_curve_derivatives(curve_geomdl, params, order=1)
    assert(allclose(geomdl_derivatives, curve.derivatives_at(params)))

    # tangents
    tangents = curve.tangents_at(params)
    rhino_tangents = [curve_rhino.TangentAt(t) for t in params]
    rhino_tangents = [[p.X, p.Y, p.Z] for p in rhino_tangents]
    geomdl_tangents = curve_tangents(geomdl_derivatives)
    assert(allclose(tangents, rhino_tangents))
    assert(allclose(tangents, geomdl_tangents))

    # curvature
    curvature = curve.curvatures_at(params)
    rhino_curvature = [curve_rhino.CurvatureAt(t) for t in params]
    rhino_curvature = [np.linalg.norm([p.X, p.Y, p.Z]) for p in rhino_curvature]
    assert(allclose(curvature, rhino_curvature))


def test_nurbs_curve():
    """Testing NurbsCurve against other libraries (geomdl, rhino3dm).
    """
    # settings
    control_points = [[816.32, 139.74, 0.0], [860.41, 238.68, 0.0], [928.17, 199.96, 0.0], [1008.83, 136.51, 0.0],
                      [1062.6, 195.66, 0.0], [1130.35, 248.36, 0.0], [1176.6, 202.11, 0.0], [1240.39, 164.93, 0.0],
                      [1315.67, 242.37, 0.0], [1396.33, 307.97, 0.0], [1466.23, 223.01, 0.0], [1520.01, 156.33, 0.0]]
    degree = 9
    params = linspace(0., 1., 10)
    weights = [0.5, 1.1, 0.7, 2., 4., 0.3, 0.1, 0.3, 0.4, 0.5, 1.4, 3.]

    # create curve
    curve = NurbsCurve(control_points, degree, weights=weights)
    curve_geomdl = create_curve(control_points, degree, curve.knot_vector, True, weights)
    curve_rhino = rhino_curve_from_curve(curve)

    # points
    points = curve.points_at(params)
    rhino_points = [curve_rhino.PointAt(t) for t in params]
    rhino_points = [[p.X, p.Y, p.Z] for p in rhino_points]
    geomdl_points = evaluate_curve(curve_geomdl, params)
    # assert(allclose(points, rhino_points))
    assert(allclose(points, geomdl_points))

    # derivatives
    geomdl_derivatives = evaluate_curve_derivatives(curve_geomdl, params, order=1)
    assert(allclose(geomdl_derivatives, curve.derivatives_at(params)))

    # tangents
    tangents = curve.tangents_at(params)
    rhino_tangents = [curve_rhino.TangentAt(t) for t in params]
    rhino_tangents = [[p.X, p.Y, p.Z] for p in rhino_tangents]
    geomdl_tangents = curve_tangents(geomdl_derivatives)
    # assert(allclose(tangents, rhino_tangents))
    assert(allclose(tangents, geomdl_tangents))

    # curvature
    curvature = curve.curvatures_at(params)
    rhino_curvature = [curve_rhino.CurvatureAt(t) for t in params]
    rhino_curvature = [np.linalg.norm([p.X, p.Y, p.Z]) for p in rhino_curvature]
    # assert(allclose(curvature, rhino_curvature))
    assert(allclose(curvature, [0.00085, 0.00076, 0.0133, 0.0069, 0.0098, 0.00037, 0.00088, 0.00337, 0.00818, 5e-05]))


if __name__ == "__main__":
    test_curve()
    test_nurbs_curve()
