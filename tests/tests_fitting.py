import os
import json
import numpy as np
from compas.geometry import allclose


from compas_nurbs.fitting_numpy import global_curve_interpolation_with_end_derivatives
from compas_nurbs.fitting_numpy import global_curve_interpolation
from compas_nurbs.knot_vectors import knot_vector_and_params

from compas_nurbs import Curve
from compas_nurbs import DATA


def read_json(filepath):
    with open(filepath, 'r') as fp:
        return json.load(fp)


def tests_rhino_compare_params_knot_vectors():
    """Compare curve interpolation with Rhino's curve interpolation.
    """
    degree = 3
    for knot_style in ['uniform', 'chord', 'chord_sqrt']:
        data = read_json(os.path.join(DATA, "rhino_curve_data_%d_%s.json" % (degree, knot_style)))
        points = np.array(data['points'])

        kvd, uk = knot_vector_and_params(points, degree, knot_style, extended=True)

        assert(allclose(kvd, data['knot_vector'], tol=0.02))

        curve = Curve(data['control_points'], data['degree'], data['knot_vector'])
        D0, Dn = curve.derivatives_at([0, 1], order=1)  # TODO: how to estimate derivatives

        control_points, knot_vector = global_curve_interpolation_with_end_derivatives(points, degree, knot_style=knot_style, periodic=False, start_derivative=D0, end_derivative=Dn)

        if knot_style != 'chord':
            assert(allclose(control_points, data['control_points'], tol=1.4))

    knot_style = 'uniform'
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


def test_end_derivatives():

    data = read_json(os.path.join(DATA, "end_derivatives_data.json"))
    points = data['points']
    start_derivative = data['start_derivative']
    end_derivative = data['end_derivative']
    knot_style = data['knot_style']
    degree = data['degree']

    control_points, knot_vector = global_curve_interpolation_with_end_derivatives(points, degree, knot_style, False, start_derivative, end_derivative)

    curve = Curve(control_points, degree, knot_vector)
    curve.to_json(os.path.join(DATA, "crv.json"))


def test_interpolation():

    data = read_json(os.path.join(DATA, "end_derivatives_data.json"))
    points = data['points']
    knot_style = data['knot_style']
    degree = data['degree']

    control_points, knot_vector = global_curve_interpolation(points, degree, knot_style, False)

    curve = Curve(control_points, degree, knot_vector)
    curve.to_json(os.path.join(DATA, "crv.json"))


if __name__ == "__main__":
    test_end_derivatives()
    test_interpolation()
    tests_rhino_compare_params_knot_vectors()
