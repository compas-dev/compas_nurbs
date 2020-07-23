import json
import compas
from compas.geometry import Primitive
from .knot_vectors import knot_vector_uniform

if not compas.IPY:
    import numpy as np
    from .evaluators_numpy import evaluate_curve
    from .evaluators_numpy import evaluate_curve_derivatives
else:
    from .evaluators import evaluate_curve
    from .evaluators import evaluate_curve_derivatives


class BSpline(Primitive):
    """
    https://github.com/orbingol/NURBS-Python/tree/5.x/geomdl

    or base on scipy's bspline?
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.BSpline.html

    BSline.basis_element(t[, extrapolate])
    BSline.derivative(self[, nu])
    BSline.integrate(self, a, b[, extrapolate])
    """

    def __init__(self, control_points, degree, knot_vector):
        self.degree = degree
        self.knot_vector = knot_vector
        if compas.IPY:
            self.control_points = control_points
        else:
            self.control_points = np.array(control_points)

    @classmethod
    def from_uniform_knot_style(cls, control_points, degree, periodic='False'):
        number_of_control_points = len(control_points)
        knot_vector = knot_vector_uniform(degree, number_of_control_points, open=True, periodic=periodic)
        return cls(control_points, degree, knot_vector)

    @property
    def is_planar(self):
        """Returns ``True`` if the curve is planar.
        """
        # TODO: only checks if the curve is planar in xy-, xz-, yz-plane, oriented bounding box to check?
        ranges = np.ptp(self.control_points, axis=0)
        is_planar = not ranges.all()
        ix, iy, iz = np.argsort(ranges)[::-1]
        return is_planar, np.take(self.control_points, [ix, iy], axis=1)

    def evaluate_at(self, params, rational=False):
        return evaluate_curve(self.control_points, self.degree, self.knot_vector, params, rational)

    def tangent_at(self, params):
        return evaluate_curve_derivatives(self.control_points, self.degree, self.knot_vector, params, order=1)

    @property
    def data(self):
        """dict: The data dictionary that represents the curve."""
        return {'control_points': [list(point) for point in self.control_points],
                'knot_vector': self.knot_vector,
                'degree': self.degree}

    @classmethod
    def from_data(cls, data):
        return cls(data['control_points'], data['degree'], data['knot_vector'])

    def to_json(self, filename):
        with open(filename, 'w') as f:
            json.dump(self.data, f, indent=4, sort_keys=True)

    @classmethod
    def from_json(cls, filename):
        with open(filename, 'r') as f:
            return cls.from_data(json.load(f))


if __name__ == "__main__":
    pass
