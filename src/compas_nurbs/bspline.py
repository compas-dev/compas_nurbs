import compas
from compas.geometry import Primitive

from compas_nurbs.knot_vectors import check_knot_vector
from compas_nurbs.knot_vectors import knot_vector_uniform
from compas_nurbs.knot_vectors import normalize_knot_vector
from compas_nurbs.utilities import prod
from compas_nurbs.utilities import reshape

if not compas.IPY:
    from collections.abc import Iterable

    import numpy as np
    from compas.geometry import transform_points_numpy
else:
    from collections import Iterable


class BSpline(Primitive):
    """A base class for rational and non-rational B-Spline geometry.

    Contains all setters and checkers.
    """

    def __init__(self, control_points, degree, knot_vector, rational, weights=None):
        self.degree = degree                                              # (degree_u, degree_v) for surfaces
        self.__rational = rational
        self.__pdim = len(degree) if isinstance(degree, Iterable) else 1
        self.control_points = control_points                              # 2d for surfaces
        self.knot_vector = knot_vector                                    # (knotvector_u, knotvector_v) for surfaces
        self.weights = weights                                            # 2d for surfaces
        self._build_backend()

    def _build_backend(self):  # needs to be overwritten by derivative classes
        raise NotImplementedError

    @property
    def rational(self):
        return self.__rational

    @property
    def domain(self):
        if self.__pdim == 1:
            return [0., 1.]
        else:
            return [[0., 1.] for _ in range(self.__pdim)]

    @property
    def control_points(self):
        return self._control_points

    @control_points.setter
    def control_points(self, control_points):
        self._control_points = control_points
        if self.__pdim == 1:
            if len(self.control_points) < self.degree + 1:
                raise ValueError("len(control_points) must be >= degree + 1")
        else:
            ok1 = len(self.degree) == len(self.count)
            ok2 = all([c >= d + 1 for c, d in zip(self.count, self.degree)])
            if not all([ok1, ok2]):
                raise ValueError("Invalid control points")

    @property
    def count(self):
        if self.__pdim == 1:
            return len(self.control_points)
        else:
            a, c = self.control_points, []  # lambda?
            for _ in range(self.__pdim):
                c.append(len(a))
                a = a[0]
            return c

    @property
    def knot_vector(self):
        """list of float : The knot vector"""
        return self._knot_vector

    @knot_vector.setter
    def knot_vector(self, knot_vector):
        if self.__pdim == 1:
            if knot_vector:
                if not check_knot_vector(knot_vector, self.count, self.degree):
                    raise ValueError("Invalid knot vector")
                self._knot_vector = normalize_knot_vector(knot_vector)
            else:
                self._knot_vector = knot_vector_uniform(self.count, self.degree)
        else:
            if knot_vector:
                ok1 = all([check_knot_vector(kv, c, d) for kv, c, d in zip(knot_vector, self.count, self.degree)])
                ok2 = len(knot_vector) == len(self.count)
                if not all([ok1, ok2]):
                    raise ValueError("Invalid knot vector")
                self._knot_vector = [normalize_knot_vector(kv) for kv in knot_vector]
            else:
                self._knot_vector = [knot_vector_uniform(c, d) for c, d in zip(self.count, self.degree)]

    @property
    def weights(self):
        return self._weights

    @weights.setter
    def weights(self, weights):
        if self.__pdim == 1:
            if weights:
                if len(weights) != self.count:
                    raise ValueError("Invalid weights! Number of weights must be equal to number of control points.")
                else:
                    self._weights = weights
            else:
                self._weights = [1. for _ in range(self.count)]
        else:
            if not weights:
                self._weights = reshape([1. for _ in range(prod(self.count))], self.count)
            else:
                self._weights = weights  # TODO check

    # ==========================================================================
    # operations
    # ==========================================================================

    def transform(self, transformation):
        xyz = np.array(self.control_points)
        xyz = transform_points_numpy(xyz.reshape(-1, xyz.shape[-1]), transformation)  # TODO more dimen?
        self.control_points = xyz.reshape(self.count + [3])
        self._build_backend()

    def get_bounding_box(self):
        raise NotImplementedError

    def trim(self):
        raise NotImplementedError

    def split(self):
        raise NotImplementedError

    def parameters_at(self):
        raise NotImplementedError

    # ==========================================================================
    # serialisation
    # ==========================================================================

    @property
    def data(self):
        """dict: The data dictionary that represents the bspline geometry."""
        return {'control_points': self.control_points,
                'degree': self.degree,
                'knot_vector': self.knot_vector,
                'rational': self.rational,
                'weights': self.weights}

    @classmethod
    def from_data(cls, data):
        return cls(data['control_points'], data['degree'], data['knot_vector'], data['rational'], data['weights'])
