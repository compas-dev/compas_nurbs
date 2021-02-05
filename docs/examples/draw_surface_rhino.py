from compas_nurbs import Surface
from compas_nurbs.rhino import SurfaceArtist

control_points = [[[0, 0, 0], [0, 4, 0.], [0, 8, -1]],
                  [[3, 0, 6], [3, 4, 0.], [3, 8, 0.]],
                  [[5, 0, 0], [5, 4, 0.], [5, 8, 3.]],
                  [[8, 0, 0], [8, 4, -3], [8, 8, 0.]]]
degree = (3, 2)

surface = Surface(control_points, degree)

artist = SurfaceArtist(surface)
artist.draw()
artist.redraw()
