import os
import json
import numpy as np
from scipy import interpolate  # noqa: F401

from compas.geometry import allclose

from compas_nurbs.fitting_numpy import global_curve_interpolation_with_end_derivatives
from compas_nurbs.fitting_numpy import global_curve_interpolation
from compas_nurbs.knot_vectors import knot_vector_and_params
from compas_nurbs.knot_vectors import CurveKnotStyle

from compas_nurbs import Curve
from compas_nurbs import DATA


def tests_rhino_compare_params_knot_vectors():
    """Compare curve interpolation with Rhino's curve interpolation.
    """

    def read_json(filepath):
        with open(filepath, 'r') as fp:
            return json.load(fp)

    degree = 3
    for knot_style in [CurveKnotStyle.Uniform, CurveKnotStyle.Chord, CurveKnotStyle.ChordSquareRoot]:
        data = read_json(os.path.join(DATA, "rhino_curve_data_%d_%s.json" % (degree, knot_style)))
        points = np.array(data['points'])

        kvd, uk = knot_vector_and_params(points, degree, knot_style, extended=True)

        assert(allclose(kvd, data['knot_vector'], tol=0.02))

        curve = Curve(data['control_points'], data['degree'], data['knot_vector'])
        D = curve.derivatives_at([0, 1], order=1)  # TODO: how to estimate derivatives
        D0, D1 = D[0, 1], D[1, 1]

        control_points, knot_vector = global_curve_interpolation_with_end_derivatives(points, degree, D0, D1, knot_style, periodic=False)

        if knot_style != CurveKnotStyle.Chord:
            assert(allclose(control_points, data['control_points'], tol=1.4))

    knot_style = CurveKnotStyle.Uniform
    for degree in [5, 7, 9]:
        data = read_json(os.path.join(DATA, "rhino_curve_data_%d_%s.json" % (degree, knot_style)))
        points = np.array(data['points'])

        kvd, uk = knot_vector_and_params(points, degree, knot_style, extended=False)

        if degree == 5:
            assert(allclose(uk, data['params'], tol=0.05))
        assert(allclose(kvd, data['knot_vector']))

        control_points, knot_vector = global_curve_interpolation(points, degree, knot_style=knot_style, periodic=False)

        if degree < 9:
            assert(allclose(control_points, data['control_points'], tol=5.8))


def test_interpolation():
    degree = 3
    points = [(0., 0., 0.), (3., 4., 0.), (-1., 4., 0.), (-4., 0., 0.), (-4., -3., 0.)]
    control_points, knot_vector = global_curve_interpolation(points, degree)
    solution = [(0, 0, 0), (6.44, 3.72, 0.), (-2.67, 7.5, 0), (-5.11, -2.72, 0), (-4, -3, 0)]
    assert(allclose(control_points, solution, tol=0.005))


def test_interpolation_with_end_derivatives():
    degree = 3
    points = [(0., 0., 0.), (3., 4., 0.), (-1., 4., 0.), (-4., 0., 0.), (-4., -3., 0.)]
    start_derivative = [17.75, 10.79, 0.0]
    end_derivative = [-0.71, -12.62, 0.0]
    result = global_curve_interpolation_with_end_derivatives(points, degree, start_derivative, end_derivative)
    solution = [(0, 0, 0), (1.48, 0.9, 0), (4.95, 5.08, 0), (-1.56, 4.87, 0), (-4.72, -0.56, 0), (-3.94, -1.95, 0), (-4, -3, 0)]
    crv = Curve(result[0], degree, result[1])
    D = crv.derivatives_at([0, 1], order=1)
    D0, D1 = D[0, 1], D[1, 1]
    assert(allclose(D0, start_derivative))
    assert(allclose(D1, end_derivative))
    assert(allclose(result[0], solution, tol=0.005))


def test_scipy_interpolation():
    filepath = os.path.join(DATA, "rhino_curve_data_3_0.json")
    with open(filepath, 'r') as fp:
        data = json.load(fp)

    points = np.array(data['points'])

    knot_style = CurveKnotStyle.Uniform
    for degree in range(2, 6): # why this runs only 2 times?
        kvd, uk = knot_vector_and_params(points, degree, knot_style, extended=False)
        control_points, knot_vector = global_curve_interpolation(points, degree, knot_style=knot_style, periodic=False)
        # why this makes such a strange behaviour for testing?
        # ...Windows fatal exception: code 0xc0000374
        # (t, c, k), u = interpolate.splprep(points.T, k=degree, task=-1, t=kvd, u=uk)
        # P = np.array(c).T
        # assert(allclose(control_points, P))


if __name__ == "__main__":
    test_scipy_interpolation()
