import compas
from compas_nurbs.bspline_curve import Curve

class NurbsCurve(Curve):

    def __init__(self, control_points, degree, knot_vector=None, weights=None):
        super(NurbsCurve, self).__init__(control_points, degree, knot_vector, weights, rational=True)


if __name__ == "__main__":
    control_points = [(0, 0, 0), (3, 4, 0), (-1, 4, 0), (-4, 0, 0), (-4, -3, 0)]
    curve = NurbsCurve(control_points, 3)
