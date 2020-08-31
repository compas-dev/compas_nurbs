============================================================
compas_nurbs: NURBS for COMPAS
============================================================

.. start-badges

.. image:: https://img.shields.io/badge/License-MIT-blue.svg
    :target: https://github.com/gramaziokohler/compas_nurbs/blob/master/LICENSE
    :alt: License MIT

.. image:: https://travis-ci.org/gramaziokohler/compas_nurbs.svg?branch=master
    :target: https://travis-ci.org/gramaziokohler/compas_nurbs
    :alt: Travis CI

.. end-badges

This package is inspirated and partly builds upon the NURBS-Python_ package, however uses a NumPy_-based backend for better performance if using CPython.
It provides wrappers about the different backends, which are streamlined with the COMPAS API.

Following Rhino3D_'s terminology, ``Curve``, and ``Surface`` are a non-uniform non-rational BSpline geometries (NUBS), ``NurbsCurve``, and ``NurbsSurface`` are 
non-uniform rational BSpline Geometries (NURBS).

Rational shapes use homogeneous coordinates which includes a weight alongside with the Cartesian coordinates.
Rational B-Splines are also named as NURBS (Non-Uniform Rational Basis Spline) and non-rational B-splines are sometimes named as NUBS (Non-Uniform Basis Spline) or directly as B-Splines.

https://github.com/nortikin/sverchok/blob/516164038e5e682b39065d657900d19dedb6ec84/utils/surface/data.py

https://github.com/nortikin/sverchok

https://github.com/alecjacobson/geometry-processing-curvature


.. _NURBS-Python: https://github.com/orbingol/NURBS-Python
.. _NumPy: https://numpy.org/
.. _Rhino3D: https://www.rhino3d.com/


Main features
-------------

* classes Curve, Surface, NurbsCurve, NurbsSurface
* Curve and Surface Interpolation
* ...

**compas_nurbs** runs on Python x.x and x.x.


Documentation
-------------

.. Explain how to access documentation: API, examples, etc.

..
.. optional sections:

Requirements
------------

.. Write requirements instructions here


Installation
------------

.. Write installation instructions here


Contributing
------------

Make sure you setup your local development environment correctly:

* Clone the `compas_nurbs <https://github.com/gramaziokohler/compas_nurbs>`_ repository.
* Install development dependencies and make the project accessible from Rhino:

::

    pip install -r requirements-dev.txt
    invoke add-to-rhino

**You're ready to start working!**

During development, use tasks on the
command line to ease recurring operations:

* ``invoke clean``: Clean all generated artifacts.
* ``invoke check``: Run various code and documentation style checks.
* ``invoke docs``: Generate documentation.
* ``invoke test``: Run all tests and checks in one swift command.
* ``invoke add-to-rhino``: Make the project accessible from Rhino.
* ``invoke``: Show available tasks.

For more details, check the `Contributor's Guide <CONTRIBUTING.rst>`_.


Releasing this project
----------------------

.. Write releasing instructions here


.. end of optional sections
..

Credits
-------------

This package was created by Romana Rust <rust@arch.ethz.ch> `@romanaru <https://github.com/romanaru>`_ at `@gramaziokohler <https://github.com/gramaziokohler>`_
