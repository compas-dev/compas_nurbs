
import json
import compas
from compas.geometry import Shape
from .knot_vectors import knot_vector_uniform
from compas_nurbs.utilities import unflatten

if not compas.IPY:
    import numpy as np
    from .evaluators_numpy import evaluate_surface
    from .evaluators_numpy import evaluate_surface_derivatives
else:
    from .evaluators import evaluate_surface
    from .evaluators import evaluate_surface_derivatives


# TODO: compas.geometry.Primitive: add from_json, to_json


class Surface(Shape):
    """
    """

    def __init__(self, control_points, count_u, count_v, degree_u, degree_v, knot_vector_u, knot_vector_v):
        if count_u * count_v != len(control_points):
            raise ValueError("u_count * v_count != len(control_points)")
        self.degree_u = degree_u
        self.degree_v = degree_v
        self.knot_vector_u = knot_vector_u
        self.knot_vector_v = knot_vector_v
        self.count_u = count_u
        self.count_v = count_v
        if compas.IPY:
            self.control_points = control_points
            self.control_points_2d = unflatten(control_points, count_u)
        else:
            self.control_points = np.array(control_points)
            self.control_points_2d = self.control_points.reshape(count_u, count_v, self.control_points[0].shape[0])

    @classmethod
    def from_uniform_knot_style(cls, control_points, count_u, count_v, degree_u, degree_v, periodic='False'):
        knot_vector_u = knot_vector_uniform(degree_u, count_u, open=True, periodic=periodic)
        knot_vector_v = knot_vector_uniform(degree_v, count_v, open=True, periodic=periodic)
        return cls(control_points, count_u, count_v, degree_u, degree_v, knot_vector_u, knot_vector_v)

    def evaluate_at(self, params_u, params_v, rational=False):
        return evaluate_surface(self.control_points_2d, self.degree_u, self.degree_v, self.knot_vector_u, self.knot_vector_v, params_u, params_v, rational)

    def tangent_at(self, params):
        return evaluate_surface_derivatives(self.control_points_2d, self.degree_u, self.degree_v, self.knot_vector_u, self.knot_vector_v, params, order=1)

    def to_obj(self):
        pass

    def transform(self, transformation):
        pass

    @property
    def data(self):
        """dict: The data dictionary that represents the surface."""
        return {'control_points': [list(point) for point in self.control_points],
                'count_u': self.count_u,
                'count_v': self.count_v,
                'degree_u': self.degree_u,
                'degree_v': self.degree_v,
                'knot_vector_u': self.knot_vector_u,
                'knot_vector_v': self.knot_vector_v, }

    @classmethod
    def from_data(cls, data):
        return cls(data['control_points'], data['count_u'], data['count_v'], data['degree_u'], data['degree_v'], data['knot_vector_u'], data['knot_vector_v'])

    def to_json(self, filename):
        with open(filename, 'w') as f:
            json.dump(self.data, f, indent=4, sort_keys=True)

    @classmethod
    def from_json(cls, filename):
        with open(filename, 'r') as f:
            return cls.from_data(json.load(f))


if __name__ == "__main__":
    pass
