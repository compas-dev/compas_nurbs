import pytest
from compas.geometry import close
from compas.geometry import allclose
from compas_nurbs import Surface


@pytest.fixture(autouse=True)
def add_imports(doctest_namespace):
    doctest_namespace["Surface"] = Surface
    doctest_namespace["close"] = close
    doctest_namespace["allclose"] = allclose


@pytest.fixture(scope='function', autouse=True)
def create_surface(request, doctest_namespace):
    if request.module.__name__ in ('compas_nurbs.bspline_surface'):
        control_points_2d = [[[0, 0, 0], [0, 4, 0.], [0, 8, -3]],
                             [[2, 0, 6], [2., 4, 0.], [2, 8, 0.]],
                             [[4, 0, 0], [4., 4, 0.], [4, 8, 3.]],
                             [[6, 0, 0], [6., 4, -3], [6, 8, 0.]]]
        degree_u, degree_v = 3, 2
        surface = Surface(control_points_2d, degree_u, degree_v)
        doctest_namespace["surface"] = surface
        yield
    else:
        yield
