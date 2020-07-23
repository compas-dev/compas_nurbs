def rhino_nurbs_surface_from_surface(surface):
    import Rhino.Geometry as rg

    point_count = surface.count_u, surface.count_v
    points = [rg.Point3d(*p) for p in surface.control_points]
    degree = surface.degree_u, surface.degree_v
    knots_u = [u for u in surface.knot_vector_u][1:-1]
    knots_v = [v for v in surface.knot_vector_v][1:-1]
    weights = None  # TODO read out

    ns = rg.NurbsSurface.Create(3, weights is not None, degree[0] + 1, degree[1] + 1, point_count[0], point_count[1])
    controlpoints = ns.Points
    index = 0
    for i in range(point_count[0]):
        for j in range(point_count[1]):
            if weights:
                cp = rg.ControlPoint(points[index], weights[index])
                controlpoints.SetControlPoint(i, j, cp)
            else:
                cp = rg.ControlPoint(points[index])
                controlpoints.SetControlPoint(i, j, cp)
            index += 1

    for i in range(ns.KnotsU.Count):
        ns.KnotsU[i] = knots_u[i]
    for i in range(ns.KnotsV.Count):
        ns.KnotsV[i] = knots_v[i]
    assert(ns.IsValid)
    return ns
