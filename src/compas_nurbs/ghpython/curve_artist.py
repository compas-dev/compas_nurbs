import Rhino.Geometry as rg

from compas_rhino.artists import BaseArtist


class CurveArtist(BaseArtist):
    """
    """

    def __init__(self, curve):
        self.curve = curve

    def draw(self):
        points = [rg.Point3d(*p) for p in self.curve.control_points]
        knots = self.curve.knot_vector[1:-1]
        cvcount = len(points)
        knotcount = cvcount + self.curve.degree - 1
        if len(knots) != knotcount:
            raise Exception("Number of elements in knots must equal the number of elements in points plus degree minus 1")
        if self.curve.weights and len(self.curve.weights) != cvcount:
            raise Exception("Number of elements in weights should equal the number of elements in points")
        rational = (self.curve.weights is not None)
        nc = rg.NurbsCurve(3, rational, self.curve.degree+1, cvcount)
        if rational:
            for i in range(cvcount):
                nc.Points.SetPoint(i, points[i], self.curve.weights[i])
        else:
            for i in range(cvcount):
                nc.Points.SetPoint(i, points[i])
        for i in range(knotcount):
            nc.Knots[i] = knots[i]
        return nc

    def draw_points(self):
        return [rg.Point3d(*p) for p in self.curve.control_points]
