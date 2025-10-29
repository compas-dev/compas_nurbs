
Changelog
=========

Unreleased
----------

**Added**

* Added ``Surface.isocurve``

**Changed**

* Migrated from ``setup.py``/``setup.cfg`` to ``pyproject.toml`` (PEP 517/518)
* Updated CI workflows to use ``compas-dev/compas-actions@v4``
* Updated tooling to use ``ruff`` for linting (replaced flake8, autopep8, isort, pydocstyle, pylint)
* Updated to ``compas_invocations2`` for invoke tasks
* Updated to ``sphinx_compas2_theme`` for documentation
* Updated repository URLs from gramaziokohler to compas-dev organization

**Fixed**

* Fixed ambiguous variable name in test_surface.py (E741)

**Deprecated**

**Removed**

* Removed deprecated files: ``setup.py``, ``setup.cfg``, ``.bumpversion.cfg``, ``pytest.ini``, ``MANIFEST.in``
* Removed ``__version__.py`` (version now defined in ``__init__.py``)

0.2.1
----------
**Changed**

* updated dependencies ``compas`` from 1.17 to 2.13

0.2.0
----------

**Added**

* Added ``SurfaceArtist`` for GhPython
* Added ``CurveArtist`` for GhPython

**Changed**

**Fixed**

**Deprecated**

**Removed**

0.1.0
-------

* Initial version
