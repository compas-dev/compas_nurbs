"""

Intro to project ...


Setup
=====

In order to use this library, ...


Main concepts
=============

Describe typical classes found in project

.. autoclass:: SampleClassName
   :members:


"""

import os
from .bspline_curve import Curve
from .bspline_surface import Surface

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "data")

__all__ = ['Curve', 'Surface']
