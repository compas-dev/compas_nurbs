import compas
from compas_nurbs import Curve
from compas_nurbs import RationalCurve
from compas_nurbs import Surface
from compas_nurbs import RationalSurface
from compas_nurbs.utilities import reshape

if compas.RHINO:
    from compas_rhino.geometry import RhinoCurve as BaseRhinoCurve
    from compas_rhino.geometry import RhinoSurface as BaseRhinoSurface

    class RhinoCurve(BaseRhinoCurve):

        @classmethod
        def from_geometry(cls, geometry):
            curve = cls()
            curve.geometry = geometry
            return curve

        def to_compas(self):
            control_points = [[p.X, p.Y, p.Z] for p in self.geometry.Points]
            knot_vector = [k for k in self.geometry.Knots]
            knot_vector = [knot_vector[0]] + knot_vector + [knot_vector[-1]]
            weights = [p.Weight for p in self.geometry.Points]
            if all([w == 1. for w in weights]):
                return Curve(control_points, self.geometry.Degree, knot_vector)
            else:
                return RationalCurve(control_points, self.geometry.Degree, knot_vector, weights=weights)

    class RhinoSurface(BaseRhinoSurface):

        @classmethod
        def from_geometry(cls, geometry):
            surface = cls()
            surface.geometry = geometry
            return surface

        def to_compas(self):
            surface = self.geometry
            knot_vector_u = list(surface.KnotsU)
            knot_vector_v = list(surface.KnotsV)
            knot_vector_u = [knot_vector_u[0]] + knot_vector_u + [knot_vector_u[-1]]
            knot_vector_v = [knot_vector_v[0]] + knot_vector_v + [knot_vector_v[-1]]
            knot_vector = knot_vector_u, knot_vector_v
            degree = surface.Degree(0), surface.Degree(1)
            points = [[p.X, p.Y, p.Z] for p in surface.Points]
            weights = [p.Weight for p in surface.Points]
            count_u = surface.Points.CountU
            count_v = surface.Points.CountV
            points = reshape(points, (count_u, count_v))
            if all([w == 1. for w in weights]):
                return Surface(points, degree, knot_vector)
            else:
                weights = reshape(weights, (count_u, count_v))
                return RationalSurface(points, degree, knot_vector, weights=weights)
