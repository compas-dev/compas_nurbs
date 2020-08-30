import compas

from compas.geometry import Primitive
from compas.geometry import Point
from compas.geometry import Vector
from compas.geometry import Frame

from compas_nurbs.knot_vectors import knot_vector_uniform
from compas_nurbs.knot_vectors import normalize_knot_vector

if not compas.IPY:
    import numpy as np
    from compas_nurbs.evaluators_numpy import create_curve
    from compas_nurbs.evaluators_numpy import evaluate_curve
    from compas_nurbs.evaluators_numpy import evaluate_curve_derivatives
    from compas_nurbs.calculations_numpy import normalize
    from compas_nurbs.fitting_numpy import interpolate_curve
else:
    from compas_nurbs.evaluators import create_curve
    from compas_nurbs.evaluators import evaluate_curve
    from compas_nurbs.evaluators import evaluate_curve_derivatives


class Curve(Primitive):
    """A base class for n-variate B-spline (non-rational) curves.

    Parameters
    ----------
    control_points : list of point
        The curve's control points.
    degree : int
        The degree of the curve.
    knot_vector : list of float, optional
        The knot vector of the curve.

    Attributes
    ----------
    data : dict
        The data representation of the curve.
    control_points : list of class:`compas.geometry.Point`
        The curve's control points.
    degree : int
        The degree of the curve.
    knot_vector : list of float
        The knot vector of the curve.

    Examples
    --------
    >>> control_points = [(0.6, 0.4, 0.0), (0.2, 2.5, 0.0), (6.0, 2.1, 0.0), (4.7, 4.5, 0.0)]
    >>> knot_vector = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    >>> curve = Curve(control_points, 3, knot_vector)

    Notes
    -----
    https://github.com/orbingol/NURBS-Python/tree/5.x/geomdl
    """

    def __init__(self, control_points, degree, knot_vector=None, weights=None, rational=False):
        self.degree = degree
        self.knot_vector = knot_vector or knot_vector_uniform(len(control_points), degree)
        self.rational = rational
        self.control_points = control_points
        self.weights = weights or [1. for _ in range(len(control_points))]
        self._curve = create_curve(self.control_points, self.degree, self.knot_vector, self.rational, self.weights)

    @property
    def knot_vector(self):
        """list of float : The knot vector of the curve."""
        return self._knot_vector

    @knot_vector.setter
    def knot_vector(self, knot_vector):
        self._knot_vector = normalize_knot_vector(knot_vector)

    @property
    def domain(self):
        return [0., 1.]

    # ==========================================================================
    # constructors
    # ==========================================================================

    @classmethod
    def from_points(cls, points, degree, knot_style=0, start_derivative=None, end_derivative=None, periodic=False):
        """Constructs an interpolated curve.

        Parameters
        ----------
        points : list of point
            The list of points on the curve we are looking for.
        degree : int
            The degree of the output parametric curve.
        start_derivative : vector, optional
            The start derivative of the curve. Defaults to ``None``.
        end_derivative : vector
            The end derivative of the curve. Defaults to ``None``.
        knotstyle : int, optional
            The knot style, either 0, 1, or 2 [uniform, chord, or chord_square_root].
            Defaults to 0, uniform.

        Returns
        -------
        :class:`Curve`
            The interpolated Curve.
        """
        cpts, kv = interpolate_curve(points, degree, knot_style, start_derivative, end_derivative, periodic)
        return cls(cpts, degree, kv)

    # ==========================================================================
    # evaluate
    # ==========================================================================

    def points_at(self, params):
        """Evaluates the curve's points at the given parametric positions.

        Parameters
        ----------
        params: list of float
            Evaluation parameters within the curve's domain of [0, 1]

        Returns
        -------
        points : list of :class:`Point`
            Point locations on curve at the parameters.

        Examples
        --------
        >>> control_points = [(0, 0, 0), (3, 4, 0), (-1, 4, 0), (-4, 0, 0), (-4, -3, 0)]
        >>> curve = Curve(control_points, 3)
        >>> curve.points_at([0.0, 0.5, 1.0])
        [Point(0.000, 0.000, 0.000), Point(-0.750, 3.000, 0.000), Point(-4.000, -3.000, 0.000)]
        """
        points = evaluate_curve(self._curve, params, self.rational)
        return [Point(*p) for p in points]

    def tangents_at(self, params):
        """Evaluates the unit tangent vector at the given parametric positions.

        Parameters
        ----------
        params: list of float

        Returns
        -------
        list of :class:`Vector`
            Unit tangent vectors of the curve at the given parametric positions.

        Examples
        --------
        >>> control_points = [(0, 0, 0), (3, 4, 0), (-1, 4, 0), (-4, 0, 0), (-4, -3, 0)]
        >>> curve = Curve(control_points, 3)
        >>> curve.tangents_at([0.0, 0.5, 1.0])
        [Vector(0.600, 0.800, 0.000), Vector(-0.868, -0.496, 0.000), Vector(0.000, -1.000, 0.000)]
        """
        D = self.derivatives_at(params, order=1)
        return [Vector(*v) for v in normalize(D[:, 1])]

    def curvatures_at(self, params):
        """Evaluates the curvature at the given parametric positions.

        Parameters
        ----------
        params: list of float

        Returns
        -------
        list of float
            Curvature values of the curve at the given parametric positions.

        Examples
        --------
        >>> from compas.geometry import allclose
        >>> control_points = [(0, 0, 0), (3, 4, 0), (-1, 4, 0), (-4, 0, 0), (-4, -3, 0)]
        >>> curve = Curve(control_points, 3)
        >>> curvature = curve.curvatures_at([0.0, 0.5])
        >>> allclose(curvature, [0.042667, 0.162835])
        True
        """
        D = self.derivatives_at(params, order=2)
        d1, d2 = D[:, 1], D[:, 2]
        k = np.linalg.norm(np.cross(d1, d2, axis=1), axis=1) / np.linalg.norm(d1, axis=1)**3  # TODO
        return list(k)

    def frames_at(self, params):
        """Evaluates the curve's frames at the given parametric positions.

        Parameters
        ----------
        params: list of float

        Returns
        -------
        list of :class:`Frame`

        Examples
        --------
        >>> control_points = [(0, 0, 0), (3, 4, 0), (-1, 4, 0), (-4, 0, 0), (-4, -3, 0)]
        >>> curve = Curve(control_points, 3)
        >>> curve.frames_at([0.5])
        [Frame(Point(-0.750, 3.000, 0.000), Vector(-0.868, -0.496, 0.000), Vector(0.496, -0.868, 0.000))]
        """
        D = self.derivatives_at(params, order=2)
        points, d1, d2 = D[:, 0], D[:, 1], D[:, 2]
        binormal = normalize(np.cross(d1, d2, axis=1))
        tangents = normalize(np.array(d1))
        normals = np.cross(binormal, tangents, axis=1)
        return [Frame(pt, xaxis, yaxis) for pt, xaxis, yaxis in zip(points, tangents, normals)]

    def osculating_circle(self, u):
        radius = 1./self.curvatures_at([u])[0]
        raise NotImplementedError

    def derivatives_at(self, params, order=1):
        return evaluate_curve_derivatives(self._curve, params, order=order, rational=self.rational)

    # ==========================================================================
    # operations
    # ==========================================================================

    def reverse(self):
        raise NotImplementedError

    def split(self):
        raise NotImplementedError

    def trim(self):
        raise NotImplementedError

    def get_bounding_box(self):
        raise NotImplementedError

    def point_at_length(self):
        raise NotImplementedError

    def parameter_at(self, points):
        raise NotImplementedError

    def transform(self, transformation):
        Point.transform_collection(self.control_points, transformation)
        self._curve = create_curve(self.control_points, self.degree, self.knot_vector, self.rational, self.weights)

    # ==========================================================================
    # queries
    # ==========================================================================

    def is_linear(self):
        raise NotImplementedError

    def is_periodic(self):
        raise NotImplementedError

    @property
    def is_planar(self):
        """Returns ``True`` if the curve is planar.
        """
        # TODO: only checks if the curve is planar in xy-, xz-, yz-plane, oriented bounding box to check?
        ranges = np.ptp(self.control_points, axis=0)
        is_planar = not ranges.all()
        ix, iy, iz = np.argsort(ranges)[::-1]
        return is_planar, np.take(self.control_points, [ix, iy], axis=1)

    # ==========================================================================
    # serialisation
    # ==========================================================================

    @property
    def data(self):
        """dict: The data dictionary that represents the curve."""
        return {'control_points': [list(point) for point in self.control_points],
                'knot_vector': self.knot_vector,
                'degree': self.degree}

    @classmethod
    def from_data(cls, data):
        return cls(data['control_points'], data['degree'], data['knot_vector'])


if __name__ == '__main__':
    import doctest
    from compas.geometry import allclose  # noqa: F401
    doctest.testmod(globs=globals())
