import numpy as np


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


def surface_normals(surface, params):
    skl = surface.derivatives_at(params, order=1)
    vectors = np.cross(skl[:, 1, 0], skl[:, 0, 1])
    return normalize_vectors(vectors)
