import rhino3dm
from compas.geometry import allclose
from compas_nurbs import BSpline
from compas_nurbs import Surface
from compas_nurbs.knot_vectors import knot_vector_uniform


def test_compare_bspline_with_rhino_curve():
    """The default rhino curve evaluates the same way as BSpline.from_uniform_knot_style.
    """
    control_points = [[816.32, 139.74, 0.0], [860.41, 238.68, 0.0], [928.17, 199.96, 0.0], [1008.83, 136.51, 0.0],
                      [1062.6, 195.66, 0.0], [1130.35, 248.36, 0.0], [1176.6, 202.11, 0.0], [1240.39, 164.93, 0.0],
                      [1315.67, 242.37, 0.0], [1396.33, 307.97, 0.0], [1466.23, 223.01, 0.0], [1520.01, 156.33, 0.0]]
    degree = 9
    n = 10
    params = [i / float(n) for i in range(n + 1)]

    # create BSpline
    curve = BSpline.from_uniform_knot_style(control_points, degree)

    # create Rhino curve
    rhino_control_points = [rhino3dm.Point3d(*p) for p in control_points]
    rhino_curve = rhino3dm.Curve.CreateControlPointCurve(rhino_control_points, degree)

    # evaluate them
    points = curve.evaluate_at(params)
    rhino_points = [rhino_curve.PointAt(t * rhino_curve.Domain.T1) for t in params]  # rhino curve domain is different
    rhino_points = [[p.X, p.Y, p.Z] for p in rhino_points]

    # compare
    assert(allclose(points, rhino_points))


def test_compare_surface_with_geomdl_surface():
    """
    """
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


if __name__ == "__main__":
    test_compare_bspline_with_rhino_curve()
