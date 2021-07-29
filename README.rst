============================================================
NURBS for COMPAS
============================================================

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

Please refer to the Documentation_ for details.

.. _NURBS-Python: https://github.com/orbingol/NURBS-Python
.. _NumPy: https://numpy.org/
.. _Documentation: https://gramaziokohler.github.io/compas_nurbs/latest/

**COMPAS NURBS** runs on Python x.x and x.x.


Getting Started
---------------

The recommended way to install **COMPAS NURBS** is to use a `Anaconda/conda <https://conda.io/docs/>`_ environment:

::

    conda config --add channels conda-forge
    conda create -n ENV_NAME compas
    conda activate ENV_NAME
    pip install compas_nurbs


Once the installation is completed, you can verify your setup.
Start Python from the command prompt and run the following:

::

    >>> import compas_nurbs

You are ready to use **COMPAS NURBS**!

Contributing
------------

Make sure you setup your local development environment correctly:

* Clone the `compas_nurbs <https://github.com/gramaziokohler/compas_nurbs>`_ repository.
* Install development dependencies and make the project accessible from Rhino:

::

    pip install -r requirements-dev.txt
    python -m compas_rhino.install

**You're ready to start working!**

During development, use tasks on the
command line to ease recurring operations:

* ``invoke clean``: Clean all generated artifacts.
* ``invoke check``: Run various code and documentation style checks.
* ``invoke docs``: Generate documentation.
* ``invoke test``: Run all tests and checks in one swift command.
* ``invoke``: Show available tasks.

For more details, check the `Contributor's Guide <CONTRIBUTING.rst>`_.


Releasing this project
----------------------

Ready to release a new version of **COMPAS NURBS**? Here's how to do it:

* We use `semver <https://semver.org/>`_, i.e. we bump versions as follows:

  * ``patch``: bugfixes.
  * ``minor``: backwards-compatible features added.
  * ``major``: backwards-incompatible changes.

* Update the ``CHANGELOG.rst`` with all novelty!
* Ready? Release everything in one command:

::

    invoke release [patch|minor|major]

* Celebrate! 💃

Credits
-------------

This package was created by Romana Rust <rust@arch.ethz.ch> `@romanarust <https://github.com/romanarust>`_ at `@gramaziokohler <https://github.com/gramaziokohler>`_
