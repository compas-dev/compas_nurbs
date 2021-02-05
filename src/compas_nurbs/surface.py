
import compas
from compas.geometry import Shape, Vector, Point

from compas_nurbs.bspline import BSpline
from compas_nurbs.curve import Curve
from compas_nurbs.curvature import SurfaceCurvature

if not compas.IPY:
    from compas_nurbs.evaluators import evaluate_surface
    from compas_nurbs.evaluators import evaluate_surface_derivatives
    from compas_nurbs.evaluators import calculate_surface_curvature
    from compas_nurbs.operations import surface_normals
    from compas_nurbs.operations import unify_curves


class Surface(BSpline, Shape):
    """A base class for n-variate B-spline (non-rational) surfaces.

    Parameters
    ----------
    control_points : list of list of point
        The surface' 2-dimensional array of control points.
    degree : tuple of int
        Surface degree for the u- and the v-direction: (u, v).
    knot_vector : tuple of list of float, optional
        Surface knot vectors. One for the u- and one for the v-direction: (ku, kv).

    Attributes
    ----------
    data : dict
        The data representation of the surface.
    control_points : list of list of point
        The surface' 2-dimensional array of control points.
    degree : tuple of int
        Surface degree for the u- and the v-direction: (u, v).
    knot_vector : tuple of list of float
        Surface knot vectors. One for the u- and one for the v-direction: (ku, kv).

    Examples
    --------
    >>> control_points_2d = [[[0, 0, 0], [0, 4, 0], [0, 8, -3]], [[2, 0, 6], [2, 4, 0], [2, 8, 0]], [[4, 0, 0], [4, 4, 0], [4, 8, 3]], [[6, 0, 0], [6, 4, -3], [6, 8, 0]]]
    >>> degree_u, degree_v = 3, 2
    >>> surface = Surface(control_points_2d, (degree_u, degree_v))
    """

    def __init__(self, control_points, degree, knot_vector=None, rational=False, weights=None):
        super(Surface, self).__init__(control_points, degree, knot_vector, rational, weights)

    def _build_backend(self):
        self._surface = self  # no numpy surface

    # ==========================================================================
    # constructors
    # ==========================================================================

    @classmethod
    def loft_from_curves(cls, curves, degree_v=3):
        """Creates a :class:`Surface` by lofting between curves.

        Parameters
        ----------
        curves: list of :class:`Curve`
            A list of curves.
        degree_v : int
            The degree of the resulting surface in v-direction.

        Returns
        -------
        :class:`Surface`
            The resulting surface.

        Examples
        --------
        >>>
        """
        curves = unify_curves(curves)
        degree_u = curves[0].degree
        knot_vector_u = curves[0].knot_vector
        degree_v = min(degree_v, len(curves) - 1)
        control_points = []
        count_u = len(curves[0].control_points)
        knot_vector_v = []
        for i in range(count_u):
            points = [crv.control_points[i] for crv in curves]
            c = Curve.from_points(points, degree_v)
            control_points.append(c.control_points)
            knot_vector_v = knot_vector_v or c.knot_vector
        # Rhino lofts into the opposite direction (u=>v)
        return cls(control_points, (degree_u, degree_v), (knot_vector_u, knot_vector_v))

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
        points = evaluate_surface(self._surface, params)
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
            A curvature object with several curvature quantities.

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
        :class:`numpy.array`
            A two-dimensional array.
        """
        return evaluate_surface_derivatives(self._surface, params, order=order)

    # ==========================================================================
    # operations
    # ==========================================================================

    # ==========================================================================
    # queries
    # ==========================================================================

    # ==========================================================================
    # serialisation
    # ==========================================================================


class RationalSurface(Surface):

    def __init__(self, control_points, degree, knot_vector=None, rational=True, weights=None):
        super(RationalSurface, self).__init__(control_points, degree, knot_vector, True, weights)

    @property
    def weighted_control_points(self):
        return [[[w * x, w * y, w * z, w] for (x, y, z), w in zip(cl, wl)] for cl, wl in zip(self.control_points, self.weights)]


if __name__ == "__main__":

    import doctest
    from compas.geometry import allclose  # noqa: F401
    from compas.geometry import close  # noqa: F401
    control_points = [[[0, 0, 0], [0., 4, 0.], [0, 8, -3]],
                      [[2, 0, 6], [2., 4, 0.], [2, 8, 0.]],
                      [[4, 0, 0], [4., 4, 0.], [4, 8, 3.]],
                      [[6, 0, 0], [6., 4, -3], [6, 8, 0.]]]
    degree = (3, 2)
    surface = Surface(control_points, degree)
    params = [(u, v) for u in [0, 0.5, 1] for v in [0, 0.5, 1]]

    doctest.testmod(globs=globals())
