from compas.geometry import normalize_vectors


def curve_tangents(derivatives):
    return normalize_vectors([d[1] for d in derivatives])
