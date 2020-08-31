"""

Intro to project ...


Setup
=====

In order to use this library, ...


Main concepts
=============

Curve
-----

Curvature

https://raw.githubusercontent.com/alecjacobson/geometry-processing-curvature/master/images/osculating-circle.gif


.. autoclass:: Curve
   :members:

Surface
-------

.. autoclass:: Surface
   :members:


NurbsCurve
----------

.. autoclass:: NurbsCurve
   :members:

NurbsSurface
------------

.. autoclass:: NurbsSurface
   :members:


"""

import os
from .bspline_curve import Curve
from .bspline_surface import Surface
from .nurbs_curve import NurbsCurve  # noqa: F401
from .nurbs_surface import NurbsSurface

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "data")

__all__ = ['Curve', 'Surface']
