import numpy as np
from .helpers import EPSILON
from .helpers import knotspan
from .knot_vectors import knot_vector_multiplicities


def normalize_vectors(vectors):
    return (vectors.T / np.linalg.norm(vectors, axis=1).reshape(1, -1)).T


def curve_tangents(derivatives):
    return normalize_vectors(derivatives[:, 1])


def curve_frames(derivatives):
    points, d1, d2 = derivatives[:, 0], derivatives[:, 1], derivatives[:, 2]
    binormal = normalize_vectors(np.cross(d1, d2, axis=1))
    tangents = normalize_vectors(np.array(d1))
    normals = np.cross(binormal, tangents, axis=1)
    return points, tangents, normals


def curve_curvatures(derivatives):
    d1, d2 = derivatives[:, 1], derivatives[:, 2]
    return np.linalg.norm(np.cross(d1, d2, axis=1), axis=1) / np.linalg.norm(d1, axis=1)**3


def curve_is_planar(curve):
    """Returns ``True`` if the curve is planar.
    """
    # TODO: only checks if the curve is planar in xy-, xz-, yz-plane, oriented bounding box to check?
    ranges = np.ptp(curve.control_points, axis=0)
    is_planar = not ranges.all()
    ix, iy, iz = np.argsort(ranges)[::-1]
    return is_planar, np.take(curve.control_points, [ix, iy], axis=1)


def unify_curves(curves):
    """
    """
    # 1. Raise degree to maximum degree
    max_degree = max([crv.degree for crv in curves])
    knot_vector = curves[0].knot_vector
    for crv in curves:
        if crv.degree != max_degree:
            raise NotImplementedError("Sorry, lofting curves of different degrees is currently not yet implemented")
        if crv.knot_vector != knot_vector:
            raise NotImplementedError("Sorry, lofting curves with different knot vectors is currently not yet implemented")
    # curve_elevate_degree(curve, final_degree)
    # unify knot vector
    return curves


def curve_elevate_degree(curve, final_degree):
    raise NotImplementedError


def surface_normals(surface, params):
    skl = surface.derivatives_at(params, order=1)
    vectors = np.cross(skl[:, 1, 0], skl[:, 0, 1])
    return normalize_vectors(vectors)


def curve_knot_refine(curve, knots2insert):
    """Insert a collection of knots on a curve.

    Corresponds to Algorithm A5.4 (Piegl & Tiller)

    Parameters
    ----------
    curve : :class:`Curve`
        The curve to insert the knots into
    knots2insert : list of float
        The knots to insert - a list of parameter positions within the curve domain.

    Returns
    -------
    :class:`Curve`
        A new curve with the knots inserted.
    """
    from compas_nurbs import Curve
    degree = curve.degree
    control_points = np.array(curve.control_points)
    knots = curve.knot_vector
    n = len(control_points) - 1
    m = n + degree + 1
    r = len(knots2insert) - 1
    a = knotspan(degree, knots2insert[0], knots)
    b = knotspan(degree, knots2insert[r], knots)
    knots_post = list()

    control_points_post = [None for _ in range(len(control_points) + len(knots2insert))]
    knots_post = [None for _ in range(len(knots) + len(knots2insert))]

    for i in range(0, a - degree + 1):
        control_points_post[i] = control_points[i]
    for i in range(b - 1, n + 1):
        control_points_post[i + r + 1] = control_points[i]
    for i in range(0, a + 1):
        knots_post[i] = knots[i]
    for i in range(b + degree, m + 1):
        knots_post[i + r + 1] = knots[i]
    i = ((b + degree) - 1)
    k = ((b + degree) + r)
    j = r
    while j >= 0:
        while knots2insert[j] <= knots[i] and i > a:
            control_points_post[k - degree - 1] = control_points[i - degree - 1]
            knots_post[k] = knots[i]
            k = (k - 1)
            i = (i - 1)
        control_points_post[k - degree - 1] = control_points_post[k - degree]
        for l in range(1, degree + 1):  # noqa E741
            ind = k - degree + l
            alfa = knots_post[k + l] - knots2insert[j]
            if abs(alfa) < EPSILON:
                control_points_post[ind - 1] = control_points_post[ind]
            else:
                alfa = alfa / (knots_post[k + l] - knots[i - degree + l])
                control_points_post[ind - 1] = alfa * control_points_post[ind - 1] + (1.0 - alfa) * control_points_post[ind]
        knots_post[k] = knots2insert[j]
        k = (k - 1)
        j = (j - 1)

    return Curve(control_points_post, degree, knots_post)


def surface_knot_refine(surface, knots2insert, direction):
    """Performs knot refinement on a Surface by inserting knots at various parameters.

    Parameters
    ----------
    surface : :class:`Surface`
        The surface to insert the knots into
    knots2insert : list of float
        The knots to insert - a list of parameter positions within the surface domain.
    direction : int
        The surface direction, either 0 (u) or 1 (v).

    Returns
    -------
    :class:`Surface`
        A new surface with the knots inserted.
    """
    from compas_nurbs import Curve
    from compas_nurbs import Surface
    degree = surface.degree[direction]
    knots = surface.knot_vector[direction]
    if direction == 0:
        control_points = surface.control_points.transpose(1, 0, 2)
    else:
        control_points = surface.control_points

    new_points = []
    for cptrow in control_points:
        c = curve_knot_refine(Curve(cptrow, degree, knots), knots2insert)
        new_points.append(c.control_points)

    new_knots = c.knot_vector
    if direction == 0:
        new_points = np.array(new_points).transpose(1, 0, 2)
        knot_vector = [new_knots, surface.knot_vector[1]]
    else:
        knot_vector = [surface.knot_vector[0], new_knots]

    return Surface(new_points, surface.degree, knot_vector)


def surface_isocurve(surface, direction, param):
    """Extract an isocurve from a surface.

    Parameters
    ----------
    surface : :class:`Surface`
        The surface.
    direction : int
        The surface direction, either 0 (u) or 1 (v).
    param : float
        The parameter at which to obtain the isocurve.

    Returns
    -------
    :class:`Curve`
        A curve in the provided direction.
    """
    from compas_nurbs import Curve
    knot_vector = surface.knot_vector[direction]
    degree = surface.degree[direction]

    surface.control_points = np.array(surface.control_points)

    knot_mults = knot_vector_multiplicities(knot_vector)
    req_knot_idx = -1
    for i in range(len(knot_mults)):
        if abs(param - knot_mults[i][0]) < EPSILON:
            req_knot_idx = i
            break
    num_knots2insert = degree + 1
    if req_knot_idx >= 0:
        num_knots2insert = (num_knots2insert - knot_mults[req_knot_idx][1])

    if num_knots2insert > 0:
        knots2insert = [param for _ in range(num_knots2insert)]
        newSrf = surface_knot_refine(surface, knots2insert, direction)
    else:
        newSrf = surface

    span = knotspan(degree, param, knot_vector)

    if abs(param - knot_vector[0]) < EPSILON:
        span = 0
    else:
        if abs(param - knot_vector[-1]) < EPSILON:
            if direction == 1:
                span = len(newSrf.controlPoints[0]) - 1
            else:
                span = len(newSrf.controlPoints) - 1

    if direction == 1:
        control_points = [row[span] for row in newSrf.control_points]
        return Curve(control_points, newSrf.degree[0], newSrf.knot_vector[0])
    else:
        return Curve(newSrf.control_points[span], newSrf.degree[1], newSrf.knot_vector[1])
