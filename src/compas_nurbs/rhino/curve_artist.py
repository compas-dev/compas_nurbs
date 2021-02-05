import time

import Rhino.Geometry as rg
import rhinoscriptsyntax as rs
import scriptcontext as sc

from compas_rhino.artists import BaseArtist


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

    def redraw(self, timeout=None):
        """Redraw the Rhino view.

        Parameters
        ----------
        timeout : float, optional
            The amount of time the artist waits before updating the Rhino view.
            The time should be specified in seconds.
            Default is ``None``.

        """
        if timeout:
            time.sleep(timeout)

        sc.doc.Views.Redraw()
