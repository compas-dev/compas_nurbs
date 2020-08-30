def rhino_nurbs_surface_from_surface(surface):
    import Rhino.Geometry as rg

    point_count = surface.count_u, surface.count_v
    #points = [rg.Point3d(*p) for p in surface.control_points]
    degree = surface.degree_u, surface.degree_v
    knots_u = [u for u in surface.knot_vector_u][1:-1]
    knots_v = [v for v in surface.knot_vector_v][1:-1]
    weights = None  # TODO read out

    ns = rg.NurbsSurface.Create(3, False, degree[0] + 1, degree[1] + 1, point_count[0], point_count[1])
    controlpoints = ns.Points
    index = 0
    for i in range(point_count[0]):
        for j in range(point_count[1]):
            p = surface.control_points_2d[i][j]
            if weights:
                cp = rg.ControlPoint(rg.Point3d(*p), weights[index])
            else:
                cp = rg.ControlPoint(rg.Point3d(*p))
            controlpoints.SetControlPoint(i, j, cp)
            index += 1

    for i in range(ns.KnotsU.Count):
        ns.KnotsU[i] = knots_u[i]
    for i in range(ns.KnotsV.Count):
        ns.KnotsV[i] = knots_v[i]
    assert(ns.IsValid)
    return ns


def rhino_curve_from_curve(curve):
    import Rhino.Geometry as rg
    P = [rg.Point3d(*p) for p in curve.control_points]
    kv = curve.knot_vector[1:-1]
    degree = curve.degree
    cvcount = len(P)
    knotcount = len(kv)
    nc = rg.NurbsCurve(3, False, degree+1, cvcount)
    for i in range(cvcount):
        nc.Points.SetPoint(i, P[i])
    for i in range(knotcount):
        nc.Knots[i] = kv[i]
    return nc
