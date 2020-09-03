
import os
from geomdl import BSpline
from geomdl import utilities
from geomdl.visualization import VisMPL
from compas_nurbs import DATA
from compas_nurbs import RationalSurface
from compas_nurbs.evaluators import create_surface


surface = RationalSurface.from_json(os.path.join(DATA, "cylinder.json"))
geomdl_surface = create_surface(surface.control_points, surface.degree, surface.knot_vector, rational=True, weights=surface.weights)


geomdl_surface.delta = 0.001

geomdl_surface.evaluate()
geomdl_surface.sample_size = 40

# Plot the control point polygon and the evaluated curve
vis_comp = VisMPL.VisSurface()
geomdl_surface.vis = vis_comp

# Don't pop up the plot window, instead save it as a PDF file
geomdl_surface.render()
