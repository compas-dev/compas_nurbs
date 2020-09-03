"""
Main concepts
=============

Curve
-----

.. autoclass:: Curve
   :members:

Surface
-------

.. autoclass:: Surface
   :members:


RationalCurve
----------

.. autoclass:: RationalCurve
   :members:

RationalSurface
------------

.. autoclass:: RationalSurface
   :members:


"""

import os
from .curve import Curve
from .curve import RationalCurve  # noqa: F401
from .surface import Surface
from .surface import RationalSurface

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "data")

__all__ = ['Curve', 'Surface', 'RationalCurve', 'RationalSurface']
