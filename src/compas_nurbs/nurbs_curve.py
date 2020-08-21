from compas_nurbs.bspline import BSpline


class NurbsCurve(BSpline):

    def __init__(self, control_points, degree, knot_vector):
        super(BSpline, self).__init__(control_points, degree, knot_vector)
        self.rational = True
