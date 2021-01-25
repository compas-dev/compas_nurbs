
from compas.geometry import Vector
from compas.geometry import Circle
from compas.geometry import Plane


class CurveCurvature(object):
    """A container class with several curve curvature quantities.
    """

    def __init__(self, curvature, frame):
        self.curvature = curvature
        self.frame = frame

    @property
    def tangent(self):
        return self.frame.xaxis

    @property
    def normal(self):
        return self.frame.yaxis

    @property
    def binormal(self):
        return self.frame.zaxis

    @property
    def radius(self):
        return 1./self.curvature

    @property
    def osculating_circle(self):
        center = self.frame.point + self.radius * self.normal
        return Circle(Plane(center, self.frame.zaxis), self.radius)


class SurfaceCurvature(object):
    """A container class with several surface curvature quantities.
    """

    def __init__(self, kappa, direction, normal, mean, gauss):
        self.normal = Vector(*normal)
        self.kappa = list(kappa)
        self.direction = [Vector(*d) for d in direction]
        self.mean = mean
        self.gauss = gauss

    @property
    def osculating_circle(self):
        raise NotImplementedError
