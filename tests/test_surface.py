import rhino3dm
from geomdl import BSpline
from geomdl import NURBS
from geomdl import compatibility

from compas.geometry import close
from compas.geometry import allclose
from compas.utilities import flatten

from compas_nurbs import Curve
from compas_nurbs import Surface
from compas_nurbs import RationalSurface
from compas_nurbs.utilities import linspace


def geomdl_surface_from_surface(surface):
    srf = NURBS.Surface() if surface.rational else BSpline.Surface()
    srf.degree_u, srf.degree_v = surface.degree
    if surface.rational:
        control_points = list(flatten(surface.control_points))
        ctrlptsw = compatibility.combine_ctrlpts_weights(control_points, list(flatten(surface.weights)))
        srf.ctrlpts_size_u, srf.ctrlpts_size_v = surface.count
        srf.ctrlptsw = ctrlptsw
    else:
        srf.ctrlpts2d = surface.control_points
    srf.knotvector_u, srf.knotvector_v = surface.knot_vector
    return srf


def rhino_surface_from_surface(surface):
    degree_u, degree_v = surface.degree
    count_u, count_v = surface.count
    knot_vector_u, knot_vector_v = surface.knot_vector
    srf_rhino = rhino3dm.NurbsSurface.Create(3, surface.rational, degree_u + 1, degree_v + 1, count_u, count_v)
    for u, l in enumerate(surface.control_points):
        for v, (x, y, z) in enumerate(l):
            if surface.rational:
                srf_rhino.Points[(u, v)] = rhino3dm.Point4d(x, y, z, surface.weights[u][v])
            else:
                srf_rhino.Points[(u, v)] = rhino3dm.Point4d(x, y, z, 1.)
    for i, v in enumerate(knot_vector_u[1:-1]):
        srf_rhino.KnotsU[i] = v
    for i, v in enumerate(knot_vector_v[1:-1]):
        srf_rhino.KnotsV[i] = v
    return srf_rhino


def test_surface():
    """Testing Surface against other libraries (geomdl, rhino3dm).
    """
    # settings
    control_points_2d = [[[0, 0, 0], [0, 4, 0.], [0, 8, -3]],
                         [[2, 0, 6], [2, 4, 0.], [2, 8, 0.]],
                         [[4, 0, 0], [4, 4, 0.], [4, 8, 3.]],
                         [[6, 0, 0], [6, 4, -3], [6, 8, 0.]]]
    degree_u = 3
    degree_v = 2
    surface = Surface(control_points_2d, (degree_u, degree_v))

    # create surfaces
    srf_geomdl = geomdl_surface_from_surface(surface)
    srf_rhino = rhino_surface_from_surface(surface)

    params_u = linspace(0., 1., 5)
    params_v = linspace(0., 1., 5)
    params = [(u, v) for u in params_u for v in params_v]

    # points
    points = surface.points_at(params)
    geomdl_points = srf_geomdl.evaluate_list(params)
    rhino_points = [srf_rhino.PointAt(*p) for p in params]
    rhino_points = [[p.X, p.Y, p.Z] for p in rhino_points]
    assert(allclose(points, geomdl_points))
    assert(allclose(points, rhino_points))

    # derivatives
    geomdl_derivatives = [srf_geomdl.derivatives(u, v, order=1) for u, v in params]
    assert(allclose(geomdl_derivatives, surface.derivatives_at(params)))

    # normals
    normals = surface.normals_at(params)
    rhino_normals = [srf_rhino.NormalAt(*p) for p in params]
    rhino_normals = [[p.X, p.Y, p.Z] for p in rhino_normals]
    geomdl_normals = [srf_geomdl.normal(p) for p in params]
    assert(allclose(normals, geomdl_normals))
    assert(allclose(normals, rhino_normals))

    # curvature
    curvature = surface.curvatures_at([(0.5, 0.5)])[0]
    rhino_curvature = {'direction': [[-0.93527, 0.3356, 0.11244], [-0.30401, -0.92441, 0.23033]], 'gauss': -
                       0.06905, 'mean': -0.12646, 'kappa': [-0.41808, 0.16516], 'normal': [0.18124, 0.18124, 0.9666]}
    assert(allclose(curvature.direction, rhino_curvature['direction'], tol=0.01))
    assert(allclose(curvature.kappa, rhino_curvature['kappa']))
    assert(allclose(curvature.normal, rhino_curvature['normal'], tol=0.01))
    assert(close(curvature.gauss, rhino_curvature['gauss']))
    assert(close(curvature.mean, rhino_curvature['mean']))


def test_rational_surface():
    """Testing NurbsSurface against other libraries.

    Unfortunately the rhino3dm Nurbs Surface does not return the correct result.
    Instead rhino_* results for comparison are copied from Rhino directly.
    """
    # settings
    control_points_2d = [[[0, 0, 0], [0, 4, 0.], [0, 8, -3]],
                         [[2, 0, 6], [2, 4, 0.], [2, 8, 0.]],
                         [[4, 0, 0], [4, 4, 0.], [4, 8, 3.]],
                         [[6, 0, 0], [6, 4, -3], [6, 8, 0.]]]
    degree = (3, 2)
    weights = [[0.2, 0.1, 0.3], [0.1, 0.7, 1.2], [1.2, 2., 0.4], [0.1, 1.1, 0.5]]
    surface = RationalSurface(control_points_2d, degree, weights=weights)

    # create surfaces
    srf_geomdl = geomdl_surface_from_surface(surface)
    # srf_rhino = rhino_surface_from_surface(surface)  # This surface is wrong

    params_u = linspace(0., 1., 5)
    params_v = linspace(0., 1., 5)
    params = [(u, v) for u in params_u for v in params_v]

    # points
    points = surface.points_at(params)
    geomdl_points = srf_geomdl.evaluate_list(params)
    rhino_points = [[0.0, 0.0, 0.0], [0.0, 1.77778, -0.33333], [0.0, 4.57143, -1.28571], [0.0, 6.85714, -2.31429], [0.0, 8.0, -3.0], [2.58947, 0.0, 0.85263],
                    [2.66072, 2.90039, 0.24467], [2.51172, 4.70621, -0.0269], [2.23804, 6.30153, -0.18823], [1.84305, 8.0, -0.30269], [3.71429, 0.0, 0.42857],
                    [3.67137, 2.70161, -0.00907], [3.52817, 4.19718, -0.07394], [3.25543, 5.67391, 0.05707], [2.67857, 8.0, 0.48214], [4.07735, 0.0, 0.14917],
                    [4.34132, 2.71191, -0.50399], [4.39635, 3.98729, -0.56116], [4.35627, 5.26156, -0.2744], [4.11864, 8.0, 0.88983], [6.0, 0.0, 0.0],
                    [6.0, 3.8, -2.475], [6.0, 4.57143, -2.35714], [6.0, 5.57143, -1.76786], [6.0, 8.0, 0.0]]
    assert(allclose(points, geomdl_points))
    assert(allclose(points, rhino_points))

    # derivatives
    derivatives = surface.derivatives_at(params, order=1)
    geomdl_derivatives = [srf_geomdl.derivatives(u, v, order=1) for u, v in params]
    assert(allclose(geomdl_derivatives, derivatives))

    # normals
    normals = surface.normals_at(params)
    geomdl_normals = [srf_geomdl.normal(p) for p in params]
    rhino_normals = [[-0.94868, 0.0, 0.31623], [-0.66692, 0.19846, 0.71822], [-0.66732, 0.27485, 0.6922], [-0.72062, 0.3195, 0.61533], [-0.76822, 0.38411, 0.51215],
                     [0.15402, 0.21472, 0.96446], [0.09377, 0.17701, 0.97973], [-0.11239, 0.1064, 0.98795], [-0.47128, -0.02486, 0.88163], [-0.77464, -0.16308, 0.61102],
                     [0.56231, 0.15519, 0.81223], [0.45653, 0.11718, 0.88196], [0.28774, 0.01626, 0.95757], [-0.05787, -0.15166, 0.98674], [-0.50649, -0.30352, 0.80706],
                     [0.53366, 0.17691, 0.82699], [0.68907, 0.05998, 0.72221], [0.63712, -0.07255, 0.76734], [0.44518, -0.25925, 0.85709], [0.00122, -0.42503, 0.90518],
                     [-0.0, 0.6, 0.8], [0.74082, 0.11858, 0.66115], [0.81509, -0.23074, 0.5314], [0.78765, -0.34549, 0.51014], [0.76822, -0.38411, 0.51215]]
    assert(allclose(normals, geomdl_normals))
    assert(allclose(geomdl_normals, rhino_normals))

    # curvature
    curvature = surface.curvatures_at([(0.5, 0.5)])[0]
    rhino_curvature = {'direction': [[-0.93542, 0.2192, 0.27737], [-0.20539, -0.97554, 0.07829]], 'gauss': -
                       0.06474, 'mean': -0.14782, 'kappa': [-0.44209, 0.14644], 'normal': [0.28774, 0.01626, 0.95757]}
    assert(allclose(curvature.direction, rhino_curvature['direction'], tol=0.01))
    assert(allclose(curvature.kappa, rhino_curvature['kappa']))
    assert(allclose(curvature.normal, rhino_curvature['normal'], tol=0.01))
    assert(close(curvature.gauss, rhino_curvature['gauss']))
    assert(close(curvature.mean, rhino_curvature['mean']))


def test_loft_surface():
    curves = []
    control_points = [(-24.265, 2.447, 0.000), (-16.799, -7.922, 0.000), (-15.554, 2.613, 0.000), (-2.198, -9.333, 0.000),
                      (-0.456, -2.447, 0.000), (13.646, 4.936, 0.000), (19.370, 0.124, 0.000), (16.550, -8.669, 0.000)]
    curves.append(Curve(control_points, 7))
    control_points = [(-22.606, -7.673, 14.030), (-18.458, 0.456, 14.030), (-8.337, -3.609, 14.030), (-5.434, 6.097, 14.030),
                      (7.673, 0.622, 14.030), (9.167, -6.512, 14.030), (19.868, -3.692, 14.030)]
    curves.append(Curve(control_points, 3))
    control_points = [(-24.845, -9.250, 28.148), (-16.052, 19.536, 28.148), (0.041, -16.799, 28.148), (3.526, 18.209, 28.148),
                      (8.503, 9.167, 28.148), (11.158, 4.687, 28.148), (18.043, -16.882, 28.148), (22.606, -2.696, 28.148)]
    curves.append(Curve(control_points, 7))
    # TODO: loft surface with curves of different degrees and knot vectors
    # surface = Surface.loft_from_curves(curves, degree_v=3)


if __name__ == "__main__":
    test_surface()
    test_rational_surface()
    test_loft_surface()
