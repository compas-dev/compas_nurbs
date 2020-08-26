"""

Intro to project ...


Setup
=====

In order to use this library, ...


Main concepts
=============

Describe typical classes found in project

.. autoclass:: Curve
   :members:


"""

import os
from .bspline_curve import Curve
from .bspline_surface import Surface
from .nurbs_curve import NurbsCurve
#from .nurbs_surface import NurbsSurface

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "data")

__all__ = ['Curve', 'Surface']
