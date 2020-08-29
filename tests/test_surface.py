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

    # settings
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

    # create Rhino surface
    srf_rhino = rhino3dm.NurbsSurface.Create(3, False, degree_u + 1, degree_v + 1, count_u, count_v)
    for u, l in enumerate(control_points_2d):
        for v, (x, y, z) in enumerate(l):
            srf_rhino.Points[(u, v)] = rhino3dm.Point4d(x, y, z, 1.)
    for i, v in enumerate(knot_vector_u[1:-1]):
        srf_rhino.KnotsU[i] = v
    for i, v in enumerate(knot_vector_v[1:-1]):
        srf_rhino.KnotsV[i] = v

    params_u = linspace(0., 1., 5)
    params_v = linspace(0., 1., 5)
    params = [(u, v) for u in params_u for v in params_v]

    # evaluate points
    geomdl_points = evaluate_surface(srf_geomdl, params)
    numpy_points = srf_numpy.points_at(params)
    rhino_points = [srf_rhino.PointAt(*p) for p in params]
    rhino_points = [[p.X, p.Y, p.Z] for p in rhino_points]

    assert(allclose(numpy_points, geomdl_points))
    assert(allclose(numpy_points, rhino_points))

    # evaluate derivatives
    res2 = evaluate_surface_derivatives(srf_geomdl, params, order=2)

    # normals
    rhino_normals = [srf_rhino.NormalAt(*p) for p in params]
    rhino_normals = [[p.X, p.Y, p.Z] for p in rhino_normals]


if __name__ == "__main__":
    test_surface()
