import numpy as np
from scipy import interpolate  # noqa: F401

from compas.geometry import allclose

from compas_nurbs import Curve
from compas_nurbs.knot_vectors import knot_vector_and_params
from compas_nurbs.knot_vectors import CurveKnotStyle


def test_interpolation():
    degree = 3
    points = [(0., 0., 0.), (3., 4., 0.), (-1., 4., 0.), (-4., 0., 0.), (-4., -3., 0.)]
    curve = Curve.from_points(points, degree, knot_style=CurveKnotStyle.Uniform)
    solution = [(0, 0, 0), (6.44, 3.72, 0.), (-2.67, 7.5, 0), (-5.11, -2.72, 0), (-4, -3, 0)]
    assert(allclose(curve.control_points, solution, tol=0.005))


def test_interpolation_with_end_derivatives():
    degree = 3
    points = [(0., 0., 0.), (3., 4., 0.), (-1., 4., 0.), (-4., 0., 0.), (-4., -3., 0.)]
    start_derivative = [17.75, 10.79, 0.0]
    end_derivative = [-0.71, -12.62, 0.0]
    curve = Curve.from_points(points, degree, knot_style=CurveKnotStyle.Uniform, start_derivative=start_derivative, end_derivative=end_derivative)
    solution = [(0, 0, 0), (1.48, 0.9, 0), (4.95, 5.08, 0), (-1.56, 4.87, 0), (-4.72, -0.56, 0), (-3.94, -1.95, 0), (-4, -3, 0)]
    D = curve.derivatives_at([0, 1], order=1)
    D0, D1 = D[0, 1], D[1, 1]
    assert(allclose(D0, start_derivative))
    assert(allclose(D1, end_derivative))
    assert(allclose(curve.control_points, solution, tol=0.005))


def tests_rhino_compare_params_knot_vectors():
    """Compare curve interpolation with Rhino's curve interpolation.
    """
    degree = 3
    points = [[-4.0, -3.0, 0.0], [-3.513204, -0.855328, 0.0], [-2.345741, 1.016886, 0.0], [-0.66136, 2.434672, 0.0], [0.97341, 1.962979, 0.0],
              [0.0, 0.0, 0.0], [0.129106, -1.755842, 0.0], [3.287345, -1.823038, 0.0], [4.916863, 0.058466, 0.0], [7.033555, 2.124761, 0.0],
              [9.267841, 2.662334, 0.0], [11.031752, 0.293654, 0.0]]

    rhino_control_points = [[-4.0, -3.0, 0.0], [-3.952, -2.268, 0.0], [-3.648, -0.776, 0.0], [-2.379, 0.994, 0.0], [-0.973, 2.821, 0.0],
                            [1.796, 2.143, 0.0], [-0.306, 0.092, 0.0], [-0.393, -2.285, 0.0], [3.758, -2.274, 0.0], [4.766, 0.167, 0.0],
                            [6.865, 2.227, 0.0], [9.458, 3.32, 0.0], [10.767, 1.242, 0.0], [11.032, 0.294, 0.0]]

    rhino_knot_vector_uniform = [0.0, 0.0, 0.0, 0.0, 0.09091, 0.18182, 0.27273, 0.36364, 0.45455, 0.54545, 0.63636, 0.72727, 0.81818, 0.90909, 1.0, 1.0, 1.0, 1.0]
    rhino_knot_vector_chord = [0.0, 0.0, 0.0, 0.0, 0.0842, 0.16868, 0.25298, 0.31813, 0.40202, 0.46943, 0.59038, 0.68568, 0.79894, 0.88692, 1.0, 1.0, 1.0, 1.0]
    rhino_knot_vector_chordsqrt = [0.0, 0.0, 0.0, 0.0, 0.08789, 0.17592, 0.26386, 0.34117, 0.4289, 0.50754, 0.61287, 0.70637, 0.80831, 0.89815, 1.0, 1.0, 1.0, 1.0]
    rhino_knot_vectors = [rhino_knot_vector_uniform, rhino_knot_vector_chord, rhino_knot_vector_chordsqrt]

    for knot_style in [CurveKnotStyle.Uniform, CurveKnotStyle.Chord, CurveKnotStyle.ChordSquareRoot]:

        rhino_knot_vector = rhino_knot_vectors[knot_style]
        kvd, uk = knot_vector_and_params(points, degree, knot_style, extended=True)
        assert(allclose(kvd, rhino_knot_vector, tol=0.02))

        curve = Curve(rhino_control_points, degree, kvd)
        D = curve.derivatives_at([0, 1], order=1)  # TODO: how to estimate derivatives for degree == 3
        D0, D1 = D[0, 1], D[1, 1]
        curve = Curve.from_points(points, degree, knot_style=knot_style, start_derivative=D0, end_derivative=D1)

        if knot_style != CurveKnotStyle.Chord:
            assert(allclose(curve.control_points, rhino_control_points, tol=1.4))


def test_scipy_interpolation():
    points = [[-4.0, -3.0, 0.0], [-3.513204, -0.855328, 0.0], [-2.345741, 1.016886, 0.0], [-0.66136, 2.434672, 0.0], [0.97341, 1.962979, 0.0],
              [0.0, 0.0, 0.0], [0.129106, -1.755842, 0.0], [3.287345, -1.823038, 0.0], [4.916863, 0.058466, 0.0], [7.033555, 2.124761, 0.0],
              [9.267841, 2.662334, 0.0], [11.031752, 0.293654, 0.0]]
    points = np.array(points)

    knot_style = CurveKnotStyle.Uniform
    for degree in range(2, 6):
        kvd, uk = knot_vector_and_params(points, degree, knot_style, extended=False)
        curve = Curve.from_points(points, degree, knot_style=knot_style)
        (t, c, k), u = interpolate.splprep(points.T, k=degree, task=-1, t=kvd, u=uk)
        P = np.array(c).T
        assert(allclose(curve.control_points, P))


if __name__ == "__main__":
    tests_rhino_compare_params_knot_vectors()
    test_interpolation()
    test_interpolation_with_end_derivatives()
    test_scipy_interpolation()
