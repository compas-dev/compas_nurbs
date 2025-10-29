"""
Main concepts
=============

Curve
-----

.. autoclass:: Curve
   :members:


RationalCurve
-------------

.. autoclass:: RationalCurve
   :members:


Surface
-------

.. autoclass:: Surface
   :members:


RationalSurface
---------------

.. autoclass:: RationalSurface
   :members:


"""

import os

from .curve import Curve
from .curve import RationalCurve        # noqa: F401
from .surface import RationalSurface
from .surface import Surface

__version__ = "0.2.0"

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "data")

__all_plugins__ = ['compas_nurbs.rhino.install']
__all__ = [
    '__version__',
    'Curve',
    'RationalCurve',
    'RationalSurface',
    'Surface'
]
