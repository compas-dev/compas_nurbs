import compas
from compas.utilities import flatten
from compas_rhino.artists import BaseArtist
from compas_nurbs.utilities import reshape
from compas_nurbs import NurbsSurface

if compas.IPY:
    import Rhino.Geometry as rg
    import rhinoscriptsyntax as rs


def surface_from_rhino_surface(rhino_surface):
    """TODO: where does this go?
    """
    S = rhino_surface
    knot_vector_u = list(S.KnotsU)
    knot_vector_v = list(S.KnotsV)
    knot_vector_u = [knot_vector_u[0]] + knot_vector_u + [knot_vector_u[-1]]
    knot_vector_v = [knot_vector_v[0]] + knot_vector_v + [knot_vector_v[-1]]
    degree = S.Degree(0), S.Degree(1)
    points = [[p.X, p.Y, p.Z] for p in S.Points]
    weights = [p.Weight for p in S.Points]
    knot_vector = [knot_vector_u, knot_vector_v]
    count_u = len(knot_vector_u) - 1 - degree[0]
    count_v = len(knot_vector_v) - 1 - degree[1]
    points = reshape(points, (count_u, count_v)) 
    weights = reshape(weights, (count_u, count_v)) 
    return NurbsSurface(points, degree, knot_vector, weights=weights)


class CurveArtist(BaseArtist):
    """
    """

    def __init__(self, curve):
        self.curve = curve

    def draw(self):
        points = [rg.Point3d(*p) for p in self.curve.control_points]
        knots = self.curve.knot_vector[1:-1]
        return rs.AddNurbsCurve(points, knots, self.curve.degree, self.curve.weights)

    def draw_points(self):
        return rs.AddPoints(self.curve.control_points)


class SurfaceArtist(BaseArtist):
    """
    """

    def __init__(self, surface):
        self.surface = surface

    def draw(self):
        point_count = self.surface.count
        points = [rg.Point3d(*p) for pl in self.surface.control_points for p in pl]
        weights = list(flatten(self.surface.weights))
        knots_u = self.surface.knot_vector[0][1:-1]
        knots_v = self.surface.knot_vector[1][1:-1]
        degree = self.surface.degree
        return rs.AddNurbsSurface(point_count, points, knots_u, knots_v, degree, weights)

    def draw_points(self):
        return rs.AddPoints(self.surface.control_points)
