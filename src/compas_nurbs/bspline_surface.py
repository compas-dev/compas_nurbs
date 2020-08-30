
import compas
from compas.geometry import Shape
from .knot_vectors import knot_vector_uniform
from compas_nurbs.utilities import unflatten

if not compas.IPY:
    import numpy as np
    from .evaluators_numpy import evaluate_surface
    from .evaluators_numpy import evaluate_surface_derivatives
else:
    from .evaluators import create_surface
    from .evaluators import evaluate_surface
    from .evaluators import evaluate_surface_derivatives


# https://mcneel.github.io/rhino3dm/python/api/Surface.html
# compute_rhino3d.Brep.CreateFromLoft(curves, start, end, loftType, closed, multiple=False)

class Surface(Shape):
    """
    """

    def __init__(self, control_points_2d, degree_u, degree_v, knot_vector_u, knot_vector_v):
        # if count_u * count_v != len(control_points):
        #    raise ValueError("u_count * v_count != len(control_points)")
        self.degree_u = degree_u
        self.degree_v = degree_v
        self.knot_vector_u = knot_vector_u
        self.knot_vector_v = knot_vector_v
        #self.count_u = count_u
        #self.count_v = count_v
        self.control_points_2d = control_points_2d
        self.rational = False

        if compas.IPY:
            #self.control_points_2d = unflatten(control_points, count_u)
            self._surface = create_surface(control_points_2d, degree_u, degree_v, knot_vector_u, knot_vector_v, rational=False, weights_u=None, weights_v=None)
        else:
            #self.control_points = np.array(control_points)
            #self.control_points_2d = self.control_points.reshape(count_u, count_v, self.control_points[0].shape[0])
            self._surface = self

    @property
    def domain(self):
        raise NotImplementedError

    @property
    def degree(self):
        return (self.degree_u, self.degree_v)

    @property
    def bounding_box(self):
        raise NotImplementedError

    # ==========================================================================
    # constructors
    # ==========================================================================

    @classmethod
    def from_points(cls, points):
        pass

    @classmethod
    def from_uniform_knot_style(cls, control_points, count_u, count_v, degree_u, degree_v):
        knot_vector_u = knot_vector_uniform(degree_u, count_u)
        knot_vector_v = knot_vector_uniform(degree_v, count_v)
        return cls(control_points, count_u, count_v, degree_u, degree_v, knot_vector_u, knot_vector_v)

    # ==========================================================================
    # evaluate
    # ==========================================================================

    def points_at(self, params):
        """Evaluates the curve's points at the given parametric positions.

        Parameters
        ----------
        params: list of tuples (u, v)
            Evaluation parameters within the curve's domain of [0, 1]

        Returns
        -------
        points : list of :class:`Point`
            Point locations on the surface at the parameters.

        Examples
        --------
        >>> 
        """
        return evaluate_surface(self._surface, params, self.rational)

    def normals_at(self, params):
        pass
        skl = obj.derivatives(uv[0], uv[1], 1)
        point = skl[0][0]
        vector = linalg.vector_cross(skl[1][0], skl[0][1])
        vector = linalg.vector_normalize(vector) if normalize else vector
        return tuple(point), tuple(vector)
        pass

    def frames_at(self, params):
        pass

    def curvatures_at(self, params):
        pass

    def isocurve_at(self,):
        raise NotImplementedError

    # ==========================================================================
    # operations
    # ==========================================================================

    def transform(self, transformation):
        pass

    def trim(self):
        raise NotImplementedError

    def split(self):
        raise NotImplementedError

    def parameters_at(self):
        raise NotImplementedError

    def derivatives_at(self, params, order=1):
        return evaluate_surface_derivatives(self._surface, params, order=order)

    # ==========================================================================
    # queries
    # ==========================================================================

    # ==========================================================================
    # serialisation
    # ==========================================================================

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

    def to_obj(self):
        pass


if __name__ == "__main__":
    pass
