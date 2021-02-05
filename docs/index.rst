****************
NURBS for COMPAS
****************

.. start-badges

.. image:: https://img.shields.io/badge/License-MIT-blue.svg
    :target: https://github.com/gramaziokohler/compas_nurbs/blob/main/LICENSE
    :alt: License MIT

.. image:: https://github.com/gramaziokohler/compas_nurbs/workflows/build/badge.svg
    :target: https://github.com/gramaziokohler/compas_nurbs/actions
    :alt: Github Actions

.. end-badges

This package is inspired by the NURBS-Python_ package, however uses a NumPy_-based backend for better performance.

``Curve``, and ``Surface`` are non-uniform non-rational B-Spline geometries (NUBS), ``RationalCurve``, and ``RationalSurface`` are
non-uniform rational B-Spline Geometries (NURBS). They all built upon the class ``BSpline``.
Coordinates have to be in 3D space (x, y, z).

.. toctree::
   :maxdepth: 3
   :titlesonly:

   Introduction <self>
   gettingstarted
   examples
   api
   contributing
   license

Indices and tables
==================

* :ref:`genindex`

.. _NURBS-Python: https://github.com/orbingol/NURBS-Python
.. _NumPy: https://numpy.org/
.. _Documentation: https://gramaziokohler.github.io/compas_nurbs/latest/
