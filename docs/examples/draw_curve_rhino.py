from compas_nurbs import Curve
from compas_nurbs.rhino import CurveArtist

control_points = [(0, 0, 0), (3, 4, 0), (-1, 4, 0),
                  (-4, 0, 0), (-4, -3, 0)]
curve = Curve(control_points, 3)

artist = CurveArtist(curve)
artist.draw()
artist.draw_points()
artist.redraw()
