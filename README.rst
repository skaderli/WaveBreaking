.. image:: https://img.shields.io/pypi/v/wavebreaking.svg
        :target: https://pypi.python.org/pypi/wavebreaking

============
WaveBreaking
============

.. image:: docs/README.gif
    :alt: wavebreaking gif
    
WaveBreaking is a python package that provides detection, classification and tracking of Rossby Wave Breaking (RWB) in weather and climate data. The detection of RWB is based on analyzing the dynamical tropopause represented by a closed contour line encircling the pole as for example in Potential Vorticity (PV) fields. By applying two different breaking indices, regions of RWB are identified and differenct characteristics of the breaking events such as area and intensity are calculated. The tracking routine provides information about the temporal evolution of the wave breaking events. Finally, the implemented plotting methods allow for a first visualization. 

The detection of RWB is based on applying a wave breaking index to the dynamical tropopause. The WaveBreaking package provides to different RWB indices:

* **Streamer Index:** The streamer index is based on work by `Wernli and Sprenger (2007)`_ (and `Sprenger et al. 2017`_). Streamers are elongated structures present on the contour line the represents the dynamical tropopause. They can be described by a pair of contour points that are close together in their geographical distance but far apart in their distance connecting the points on the contour. Further description can be found `here <https://occrdata.unibe.ch/students/theses/msc/406.pdf>`_.

* **Overturning Index:** The overturning index is based on work by `Barnes and Hartmann (2012)`_. This index indetifies overturning structures of the contour line that represents the dynamical tropopause. An overturning of the contour line is present if the contour intersects at least three times with the same longitude. Further description can be found `here <https://occrdata.unibe.ch/students/theses/msc/406.pdf>`_.

.. _`Wernli and Sprenger (2007)`: https://journals.ametsoc.org/view/journals/atsc/64/5/jas3912.1.xml
.. _`Sprenger et al. 2017`: https://journals.ametsoc.org/view/journals/bams/98/8/bams-d-15-00299.1.xml
.. _`Barnes and Hartmann (2012)`: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2012JD017469

The tool is designed to analyze gridded data provided as an `xarray.Dataset <https://docs.xarray.dev/en/stable/generated/xarray.Dataset.html>`_. Output is provided either in a `pandas.DataFrame <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html>`_ or in an `xarray.Dataset <https://docs.xarray.dev/en/stable/generated/xarray.Dataset.html>`_.

The class structure and parts of the tracking routines are based on the `ConTrack - Contour Tracking <https://github.com/steidani/ConTrack>`_ tool developed by `Daniel Steinfeld <https://github.com/steidani>`_. 

Important information:
-----------------

* The package is still under constructuion and therefore, major errors could occur. 
* Free software: MIT license
* Documentation: https://wavebreaking.readthedocs.io.
* Further documentation about the implemented methods can be found `here <https://occrdata.unibe.ch/students/theses/msc/406.pdf>`_

Installation
--------
--> WILL BE ADDED SOON

Tutorial
--------
--> WILL BE ADDED SOON  

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
