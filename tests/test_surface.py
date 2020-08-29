import compas
import numpy as np
from compas.geometry._core.basic import allclose
from compas_nurbs.bspline_surface import Surface
from compas_nurbs.evaluators import create_surface, evaluate_surface, evaluate_surface_derivatives
from compas.utilities import flatten
from compas_nurbs.utilities import unflatten
from compas_nurbs.utilities import linspace
from compas_nurbs import Surface
import rhino3dm


def test_surface():

    control_points_2d = [[[0, 0, 0], [0, 4, 0.], [0, 8, -3]],
                         [[2, 0, 6], [2, 4, 0.], [2, 8, 0.]],
                         [[4, 0, 0], [4, 4, 0.], [4, 8, 3.]],
                         [[6, 0, 0], [6, 4, -3], [6, 8, 0.]]]

    degree_u = 3
    degree_v = 2

    knot_vector_u = [0, 0, 0, 0, 1, 1, 1, 1]
    knot_vector_v = [0, 0, 0, 1, 1, 1]

    srf_geomdl = create_surface(control_points_2d, degree_u, degree_v, knot_vector_u, knot_vector_v)
    srf_numpy = Surface(control_points_2d, degree_u, degree_v, knot_vector_u, knot_vector_v)

    count_v, count_u = len(control_points_2d[0]), len(control_points_2d)

    print(np.array(control_points_2d).shape)
    control_points = list(flatten(control_points_2d))
    control_points = np.array(control_points)
    control_points_2d = control_points.reshape(count_u, count_v, control_points[0].shape[0])

    params_u = linspace(0., 1., 5)
    params_v = linspace(0., 1., 5)
    params = [(u, v) for u in params_u for v in params_v]

    res2 = evaluate_surface(srf_geomdl, params)
    res1 = srf_numpy.points_at(params)

    print(res2)

    assert(allclose(res1, res2))

    res2 = evaluate_surface_derivatives(srf_geomdl, params, order=2)
    print(np.array(res2).shape)

    ns = rhino3dm.NurbsSurface.Create(3, False, degree_u, degree_v, count_u, count_v)
    P = [rhino3dm.Point3d(*p) for p in control_points]
    kv_u = knot_vector_u[1:-1]
    kv_v = knot_vector_v[1:-1]

    print(count_u, count_v)

    for u, l in enumerate(control_points_2d):
        for v, (x, y, z) in enumerate(l):
            ns.Points[(u, v)] = rhino3dm.Point4d(x, y, z, 1.)

    print(ns.Points.CountU)
    print(ns.Points.CountV)
    print(">>", len(ns.KnotsU))
    print(">>", len(ns.KnotsV))

    """
    for i in range(len(kv_u)):
        print("i", i)
        ns.KnotsU[i] = kv_u[i]
    for i in range(len(kv_v)):
        ns.KnotsV[i] = kv_v[i]
    print(dir(ns))
    """


if __name__ == "__main__":
    test_surface()
