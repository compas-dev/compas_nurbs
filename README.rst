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

This package is inspirated by the NURBS-Python_ package, however uses a NumPy_-based backend for better performance.

``Curve``, and ``Surface`` are non-uniform non-rational B-Spline geometries (NUBS), ``RationalCurve``, and ``RationalSurface`` are 
non-uniform rational B-Spline Geometries (NURBS). They all built upon the class ``BSpline``.
Coordinates have to be in 3D space (x, y, z).

Please refer to the Documentation_ for details.

.. _NURBS-Python: https://github.com/orbingol/NURBS-Python
.. _NumPy: https://numpy.org/
.. _Documentation: https://gramaziokohler.github.io/compas_nurbs/latest/

**compas_nurbs** runs on Python x.x and x.x.


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

This package was created by Romana Rust <rust@arch.ethz.ch> `@romanarust <https://github.com/romanarust>`_ at `@gramaziokohler <https://github.com/gramaziokohler>`_
