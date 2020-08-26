from compas_nurbs.bspline_curve import Curve


class NurbsCurve(Curve):

    def __init__(self, control_points, degree, knot_vector):
        super(Curve, self).__init__(control_points, degree, knot_vector)
        self.rational = True
