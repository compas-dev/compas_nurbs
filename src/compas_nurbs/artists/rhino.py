import compas
from compas.utilities import flatten
from compas_rhino.artists import BaseArtist

if compas.IPY:
    import Rhino.Geometry as rg
    import rhinoscriptsyntax as rs

class CurveArtist(BaseArtist):
    """
    """
    def __init__(self, curve):
        self.curve = curve

    def draw(self):
        points = [rg.Point3d(*p) for p in self.curve.control_points]
        knots = self.curve.knot_vector[1:-1]
        return rs.AddNurbsCurve(points, knots, self.curve.degree, self.curve.weights)

class SurfaceArtist(BaseArtist):
    """
    """
    def __init__(self, surface):
        self.surface = surface
    
    def draw(self):
        point_count = self.surface.count_u, self.surface.count_v
        points = [rg.Point3d(*p) for pl in self.surface.control_points for p in pl]
        weights = flatten(surface.weights) if surface.weights else None
        knots_u = self.surface.knot_vector_u[1:-1]
        knots_v = self.surface.knot_vector_v[1:-1]
        degree = [self.surface.degree_u, self.surface.degree_v]
        return rs.AddNurbsSurface(point_count, points, knots_u, knots_v, degree, weights)

    
