import os
import compas
import numpy as np
from compas.geometry import allclose
from compas_nurbs.bspline_surface import Surface
from compas_nurbs.evaluators import create_surface, evaluate_surface, evaluate_surface_derivatives
from compas.utilities import flatten
from compas_nurbs.utilities import unflatten
from compas_nurbs.utilities import linspace
from compas_nurbs.knot_vectors import knot_vector_uniform
from compas_nurbs import Surface
import rhino3dm
from geomdl.knotvector import generate
from compas_nurbs import DATA


def test_surface():
    # settings
    control_points_2d = [[[0, 0, 0], [0, 4, 0.], [0, 8, -3]],
                         [[2, 0, 6], [2, 4, 0.], [2, 8, 0.]],
                         [[4, 0, 0], [4, 4, 0.], [4, 8, 3.]],
                         [[6, 0, 0], [6, 4, -3], [6, 8, 0.]]]
    degree_u = 3
    degree_v = 2
    surface = Surface(control_points_2d, degree_u, degree_v)
    knot_vector_u = surface.knot_vector_u
    knot_vector_v = surface.knot_vector_v
    count_u, count_v = surface.count_u, surface.count_v

    surface.to_json(os.path.join(DATA, "surface.json"))

    # create surfaces

    srf_geomdl = create_surface(control_points_2d, degree_u, degree_v, knot_vector_u, knot_vector_v)
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
    points = surface.points_at(params)
    geomdl_points = evaluate_surface(srf_geomdl, params)
    rhino_points = [srf_rhino.PointAt(*p) for p in params]
    rhino_points = [[p.X, p.Y, p.Z] for p in rhino_points]

    assert(allclose(points, geomdl_points))
    assert(allclose(points, rhino_points))

    # evaluate derivatives
    geomdl_derivatives = evaluate_surface_derivatives(srf_geomdl, params, order=1)

    assert(allclose(geomdl_derivatives, surface.derivatives_at(params)))

    # normals
    normals = surface.normals_at(params)
    rhino_normals = [srf_rhino.NormalAt(*p) for p in params]
    rhino_normals = [[p.X, p.Y, p.Z] for p in rhino_normals]
    geomdl_normals = [srf_geomdl.normal(p)[1] for p in params]

    assert(allclose(normals, geomdl_normals))
    assert(allclose(normals, rhino_normals))

    print(surface.curvatures_at([(0.5, 0.5)]))


if __name__ == "__main__":
    test_surface()
