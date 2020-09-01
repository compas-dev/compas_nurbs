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
from .curve import Curve
from .curve import NurbsCurve  # noqa: F401
from .surface import Surface
from .surface import NurbsSurface

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "data")

__all__ = ['Curve', 'Surface', 'NurbsCurve', 'NurbsSurface']
