
import compas
from compas.geometry import Shape, Vector, Point
from geomdl.knotvector import generate

if not compas.IPY:
    from compas_nurbs.evaluators_numpy import evaluate_surface
    from compas_nurbs.evaluators_numpy import evaluate_surface_derivatives
    from compas_nurbs.evaluators_numpy import calculate_surface_curvature
    from compas_nurbs.operations_numpy import surface_normals
else:
    from compas_nurbs.evaluators import create_surface
    from compas_nurbs.evaluators import evaluate_surface
    from compas_nurbs.evaluators import evaluate_surface_derivatives


# https://mcneel.github.io/rhino3dm/python/api/Surface.html


class SurfaceCurvature(object):
    def __init__(self, kappa, direction, normal, mean, gauss):
        self.normal = Vector(*normal)
        self.kappa = list(kappa)
        self.direction = [Vector(*d) for d in direction]
        self.mean = mean
        self.gauss = gauss
        self.osculating_circle = None


class Surface(Shape):
    """A base class for n-variate B-spline (non-rational) surfaces.

    Parameters
    ----------
    control_points_2d : list of list of point
        The surface' 2-dimensional array of control points.
    degree_u : int
        Surface degree in u-direction.
    degree_v : int
        Surface degree in v-direction.
    knot_vector_u : list of float, optional
        Surface knot vector for the u-direction.
    knot_vector_v : list of float, optional
        Surface knot vector for the v-direction.

    Attributes
    ----------
    data : dict
        The data representation of the surface.
    control_points_2d : list of list of point
        The surface' 2-dimensional array of control points.
    degree_u : int
        Surface degree in u-direction.
    degree_v : int
        Surface degree in v-direction.
    knot_vector_u : list of float, optional
        Surface knot vector for the u-direction.
    knot_vector_v : list of float, optional
        Surface knot vector for the v-direction.

    Examples
    --------
    >>> control_points_2d = [[[0, 0, 0], [0, 4, 0], [0, 8, -3]], [[2, 0, 6], [2, 4, 0], [2, 8, 0]], [[4, 0, 0], [4, 4, 0], [4, 8, 3]], [[6, 0, 0], [6, 4, -3], [6, 8, 0]]]
    >>> degree_u, degree_v = 3, 2
    >>> surface = Surface(control_points_2d, degree_u, degree_v)
    """

    def __init__(self, control_points_2d, degree_u, degree_v, knot_vector_u=None, knot_vector_v=None):
        self.control_points_2d = control_points_2d
        self.degree_u = degree_u
        self.degree_v = degree_v
        # TODO checkers and setters
        self.knot_vector_u = knot_vector_u or generate(degree_u, self.count_u)
        self.knot_vector_v = knot_vector_v or generate(degree_v, self.count_v)
        self.rational = False
        if compas.IPY:
            self._surface = create_surface(control_points_2d, degree_u, degree_v, knot_vector_u, knot_vector_v, rational=False, weights_u=None, weights_v=None)
        else:
            self._surface = self  # no numpy surface

    @property
    def count_u(self):
        """The number of control points in the u-direction."""
        return len(self.control_points_2d)

    @property
    def count_v(self):
        """The number of control points in the v-direction."""
        return len(self.control_points_2d[0])

    @property
    def domain(self):
        raise NotImplementedError

    @property
    def bounding_box(self):
        raise NotImplementedError

    # ==========================================================================
    # constructors
    # ==========================================================================
    @classmethod
    def from_curves(cls, curves):
        # compute_rhino3d.Brep.CreateFromLoft(curves, start, end, loftType, closed, multiple=False)
        pass

    @classmethod
    def from_points(cls, points):
        pass

    # ==========================================================================
    # evaluate
    # ==========================================================================

    def points_at(self, params):
        """Evaluates the surface's points at the given parametric positions.

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
        >>> params = [(0.1, 0.1), (0.1, 0.5), (0.5, 0.1), (0.5, 0.5)]
        >>> surface.points_at(params)
        [Point(0.600, 0.800, 1.159), Point(0.600, 4.000, -0.164), Point(3.000, 0.800, 1.763), Point(3.000, 4.000, 0.562)]
        """
        points = evaluate_surface(self._surface, params, self.rational)
        return [Point(*p) for p in points]

    def normals_at(self, params):
        """Evaluates the surface's normals at the given parametric positions.

        Parameters
        ----------
        params: list of tuples (u, v)
            Evaluation parameters within the curve's domain of [0, 1]

        Returns
        -------
        vector : list of :class:`Vector`
            Normalized surface normals at the parameters.

        Examples
        --------
        >>> params = [(0.1, 0.1), (0.1, 0.5), (0.5, 0.1), (0.5, 0.5)]
        >>> surface.normals_at(params)
        [Vector(-0.822, 0.203, 0.533), Vector(-0.605, 0.324, 0.727), Vector(0.503, 0.424, 0.753), Vector(0.181, 0.181, 0.967)]
        """
        normals = surface_normals(self._surface, params)
        return [Vector(*n) for n in normals]

    def curvatures_at(self, params):
        """Evaluates the surface' curvature at the given parametric positions.

        Parameters
        ----------
        params: list of tuples (u, v)
            Evaluation parameters within the curve's domain of [0, 1]

        Returns
        -------
        :class:`SurfaceCurvature`
            A container class with several surface curvature quantities, such as 
            - principal curvature directions
            - principal curvature values (kappa)
            - mean curvature
            - gauss curvature
            - normal

        Examples
        --------
        >>> params = [(0.5, 0.5)]
        >>> curvature = surface.curvatures_at(params)[0]
        >>> curvature.direction # principal curvature directions
        [Vector(-0.935, 0.336, 0.112), Vector(-0.304, -0.924, 0.230)]
        >>> allclose(curvature.kappa, [-0.41808, 0.16516]) # principal curvature values
        True
        >>> curvature.normal
        Vector(0.181, 0.181, 0.967)
        >>> close(curvature.gauss, -0.06905)
        True
        >>> close(curvature.mean, -0.12646)
        True
        """
        derivatives = self._surface.derivatives_at(params, order=2)
        kappa1, kappa2, direction1, direction2, normal, mean, gauss = calculate_surface_curvature(derivatives)
        return [SurfaceCurvature((k1, k2), (d1, d2), n, m, g) for k1, k2, d1, d2, n, m, g in zip(kappa1, kappa2, direction1, direction2, normal, mean, gauss)]

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
        """Evaluates n-th order surface derivatives at the given (u, v) parameter pairs.

        Parameters
        ----------
        params : list of (u, v) tuples
            The parameters to evaluate in the [0, 1] domain.
        order : int
            The derivative order.

        Returns
        -------
        list of list of list
        """
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
        return {'control_points_2d': self.control_points_2d,
                'degree_u': self.degree_u,
                'degree_v': self.degree_v,
                'knot_vector_u': self.knot_vector_u,
                'knot_vector_v': self.knot_vector_v, }

    @classmethod
    def from_data(cls, data):
        return cls(data['control_points_2d'], data['degree_u'], data['degree_v'], data['knot_vector_u'], data['knot_vector_v'])

    def to_obj(self):
        raise NotImplementedError


if __name__ == "__main__":

    import doctest
    from compas.geometry import allclose  # noqa: F401
    from compas.geometry import close  # noqa: F401
    control_points_2d = [[[0, 0, 0], [0, 4, 0.], [0, 8, -3]],
                         [[2, 0, 6], [2., 4, 0.], [2, 8, 0.]],
                         [[4, 0, 0], [4., 4, 0.], [4, 8, 3.]],
                         [[6, 0, 0], [6., 4, -3], [6, 8, 0.]]]
    degree_u, degree_v = 3, 2
    surface = Surface(control_points_2d, degree_u, degree_v)
    params = [(u, v) for u in [0, 0.5, 1] for v in [0, 0.5, 1]]

    doctest.testmod(globs=globals())
