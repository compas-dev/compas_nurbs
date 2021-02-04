********************************************************************************
Getting Started
********************************************************************************

.. toctree::
   :maxdepth: 1
   :titlesonly:
   :glob:

Installation
============

On Mac
------

.. code-block:: bash

    conda config --add channels conda-forge
    conda create -n ENV_NAME compas python.app
    conda activate ENV_NAME
    pip install compas_nurbs

On Windows
----------

.. code-block:: bash

    conda config --add channels conda-forge
    conda create -n ENV_NAME compas
    conda activate ENV_NAME
    pip install compas_nurbs


Once the installation is completed, you can verify your setup.
Start Python from the command prompt and run the following:

::

    >>> import compas_nurbs

You are ready to use **COMPAS NURBS**!
