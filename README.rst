.. image:: https://img.shields.io/pypi/v/wavebreaking.svg
        :target: https://pypi.python.org/pypi/wavebreaking

.. image:: https://img.shields.io/travis/skaderli/wavebreaking.svg
        :target: https://travis-ci.com/skaderli/wavebreaking

.. image:: https://readthedocs.org/projects/wavebreaking/badge/?version=latest
        :target: https://wavebreaking.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

============
WaveBreaking
============

WaveBreaking is a python package that provides detection, classification and tracking of Rossby Wave Breaking (RWB) in weather and climate data. The detection of RWB is based on analyzing the dynamical tropopause represented by a closed contour line encircling the pole. By applying two different breaking indeces, regions of RWB are identified and differenct characteristics of the events such as size and intensity are calculated. Finally, the tracking method provides information about the temporal evolution of the wave breaking events.

The detection of RWB is based on two different wave breaking indices:

* Streamer index by Wernli and Sprenger (2007)
* Overturning index by Barnes and Hartmann (2012)

The tool is based on xarray...

* Free software: MIT license
* Documentation: https://wavebreaking.readthedocs.io.


Features
--------

* TODO

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
