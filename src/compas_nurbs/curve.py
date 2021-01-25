import compas

from compas.geometry import Point
from compas.geometry import Vector
from compas.geometry import Frame

from compas_nurbs.bspline import BSpline
from compas_nurbs.curvature import CurveCurvature

if not compas.IPY:
    from compas_nurbs.evaluators import create_curve
    from compas_nurbs.evaluators import evaluate_curve
    from compas_nurbs.evaluators import evaluate_curve_derivatives
    from compas_nurbs.operations import curve_tangents
    from compas_nurbs.operations import curve_frames
    from compas_nurbs.operations import curve_curvatures
    from compas_nurbs.fitting import interpolate_curve


class Curve(BSpline):
    """A non-uniform non-rational B-Spline curve.

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
    >>> degree = 3
    >>> curve = Curve(control_points, degree, knot_vector)
    """

    def __init__(self, control_points, degree, knot_vector=None, rational=False, weights=None):
        super(Curve, self).__init__(control_points, degree, knot_vector, rational, weights)

    def _build_backend(self):
        if not compas.IPY:
            self._curve = create_curve(self.control_points, self.degree, self.knot_vector, self.rational, self.weights)

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
            The interpolated curve.

        Examples
        --------
        >>> points = [(0, 0, 0), (3, 4, 0), (-1, 4, 0), (-4, 0, 0), (-4, -3, 0)]
        >>> curve = Curve.from_points(points, 3, knot_style=0)

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
        >>> curve.points_at([0.0, 0.5, 1.0])
        [Point(0.000, 0.000, 0.000), Point(-0.750, 3.000, 0.000), Point(-4.000, -3.000, 0.000)]
        """
        points = evaluate_curve(self, params)
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
        >>> curve.tangents_at([0.0, 0.5, 1.0])
        [Vector(0.600, 0.800, 0.000), Vector(-0.868, -0.496, 0.000), Vector(0.000, -1.000, 0.000)]
        """
        derivatives = self.derivatives_at(params, order=1)
        tangents = curve_tangents(derivatives)
        return [Vector(*v) for v in tangents]

    def curvatures_at(self, params):
        """Evaluates the curvature at the given parametric positions.

        Parameters
        ----------
        params: list of float

        Returns
        -------
        :class:`CurveCurvature`
            A curvature object with several curvature quantities.

        Examples
        --------
        >>> curvature = curve.curvatures_at([0.5])[0]
        >>> close(curvature.curvature, 0.162835)
        True
        >>> curvature.osculating_circle
        Circle(Plane(Point(2.297, -2.332, 0.000), Vector(0.000, 0.000, 1.000)), 6.141172894211785)
        """
        derivatives = self.derivatives_at(params, order=2)
        curvatures = curve_curvatures(derivatives)
        frames = self.frames_at(params)
        return [CurveCurvature(curvature, frame) for curvature, frame in zip(curvatures, frames)]

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
        >>> curve.frames_at([0.5])
        [Frame(Point(-0.750, 3.000, 0.000), Vector(-0.868, -0.496, 0.000), Vector(0.496, -0.868, 0.000))]
        """
        derivatives = self.derivatives_at(params, order=2)
        points, tangents, normals = curve_frames(derivatives)
        return [Frame(pt, xaxis, yaxis) for pt, xaxis, yaxis in zip(points, tangents, normals)]

    def derivatives_at(self, params, order=1):
        """Evaluates the n-th order curve derivatives at the given parametric positions.

        The output of this method is a list of n-th order derivatives. If order is 0,
        then it will only output the evaluated point. Similarly, if order is 2, then
        it will output the evaluated point, the 1st derivative and the 2nd derivative.

        Parameters
        ----------
        params: list of float
        order: int
            The derivative order.

        Returns
        -------
        :class:`numpy.array`
            A two-dimensional array.

        Examples
        --------
        >>> curve.derivatives_at([0.5], order=0)
        array([[[-0.75,  3.  ,  0.  ]]])
        >>> curve.derivatives_at([0.5], order=1)
        array([[[ -0.75,   3.  ,   0.  ],
                [-10.5 ,  -6.  ,   0.  ]]])
        >>> curve.derivatives_at([0.5], order=2)
        array([[[ -0.75,   3.  ,   0.  ],
                [-10.5 ,  -6.  ,   0.  ],
                [  6.  , -24.  ,   0.  ]]])
        """
        return evaluate_curve_derivatives(self, params, order=order)

    # ==========================================================================
    # operations
    # ==========================================================================

    def reverse(self):
        """Reverses the curve.

        Returns
        -------
        None

        Examples
        --------
        >>> control_points = [(0, 0, 0), (3, 4, 0), (-1, 4, 0), (-4, 0, 0), (-4, -3, 0)]
        >>> curve = Curve(control_points, 3)
        >>> curve.reverse()
        """
        self.control_points = list(reversed(self.control_points))
        self.weights = list(reversed(self.weights))
        self._build_backend()

    # ==========================================================================
    # queries
    # ==========================================================================


class RationalCurve(Curve):
    """A non-uniform rational B-Spline curve (NURBS-Curve).


    Examples
    --------
    >>> control_points = [(0, 0, 0), (3, 4, 0), (-1, 4, 0), (-4, 0, 0), (-4, -3, 0)]
    >>> weights = [0.3, 0.2, 1., 0.4, 2.]
    >>> curve = Curve(control_points, 3, weights=weights)
    """

    def __init__(self, control_points, degree, knot_vector=None, weights=None):
        super(RationalCurve, self).__init__(control_points, degree, knot_vector, rational=True, weights=weights)

    @property
    def weighted_control_points(self):
        """The weighted control points."""
        return [[w * x, w * y, w * z, w] for (x, y, z), w in zip(self.control_points, self.weights)]


if __name__ == '__main__':
    import doctest
    from compas.geometry import allclose  # noqa: F401
    from compas.geometry import close  # noqa: F401
    control_points = [(0, 0, 0), (3, 4, 0), (-1, 4, 0), (-4, 0, 0), (-4, -3, 0)]
    curve = Curve(control_points, 3)
    doctest.testmod(globs=globals())

    control_points = [(0, 0, 0), (3, 4, 0), (-1, 4, 0), (-4, 0, 0), (-4, -3, 0)]
    curve = RationalCurve(control_points, 3)
