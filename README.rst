.. image:: https://img.shields.io/pypi/v/wavebreaking.svg
        :target: https://pypi.python.org/pypi/wavebreaking
        
.. image:: https://img.shields.io/github/license/skaderli/wavebreaking
        :target: https://github.com/skaderli/wavebreaking/blob/master/LICENSE
        :alt: License
        
.. image:: https://readthedocs.org/projects/wavebreaking/badge/?version=latest
        :target: https://wavebreaking.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status
        
.. image:: https://www.codefactor.io/repository/github/skaderli/wavebreaking/badge
   :target: https://www.codefactor.io/repository/github/skaderli/wavebreaking
   :alt: CodeFactor

====================================================================================
WaveBreaking - Detection, Classification and Tracking of Rossby Wave Breaking
====================================================================================

.. image:: https://raw.githubusercontent.com/skaderli/WaveBreaking/main/docs/figures/readme.gif
    :alt: readme gif

.. start_intro
        
WaveBreaking is a Python package that provides detection, classification and tracking of Rossby Wave Breaking (RWB) in weather and climate data. The detection of RWB is based on analyzing the dynamical tropopause represented by a closed contour line encircling the pole as for example the 2 Potential Vorticity Units (PVU) contour line in Potential Vorticity (PV) fields. By applying three different breaking indices, regions of RWB are identified and different characteristics of the events such as area and intensity are calculated. The event tracking provides information about the temporal evolution of the RWB events. Finally, the implemented plotting methods allow for a first visualization. This tool was developed during my master studies at the University of Bern. 

The detection of RWB is based on applying a RWB index to the dynamical tropopause. The WaveBreaking package provides three different RWB indices:

* **Streamer Index:** The streamer index is based on work by `Wernli and Sprenger (2007)`_ (and `Sprenger et al. 2017`_). Streamers are elongated structures present on the contour line that represents the dynamical tropopause. They can be described by a pair of contour points that are close together considering their geographical distance but far apart considering their distance connecting the points on the contour. Further description can be found in my `master thesis <https://occrdata.unibe.ch/students/theses/msc/406.pdf>`_.

* **Overturning Index:** The overturning index is based on work by `Barnes and Hartmann (2012)`_. This index identifies overturning structures of the contour line. An overturning of the contour line is present if the contour intersects at least three times with the same longitude. Further description can be found in my `master thesis <https://occrdata.unibe.ch/students/theses/msc/406.pdf>`_.

* **Cutoff Index:** The Cutoff Index provides information about the decaying of a wave breaking event. From a Potential Vorticity perspective, a wave breaking event is formed by an elongation of the 2 PVU contour line. These so-called streamers can elongate further until they separate from the main stratospheric or tropospheric body. The separated structure is referred to as a cutoff (`Wernli and Sprenger (2007)`_.

.. _`Wernli and Sprenger (2007)`: https://journals.ametsoc.org/view/journals/atsc/64/5/jas3912.1.xml
.. _`Sprenger et al. 2017`: https://journals.ametsoc.org/view/journals/bams/98/8/bams-d-15-00299.1.xml
.. _`Barnes and Hartmann (2012)`: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2012JD017469

The tool is designed to analyze gridded data provided as an `xarray.DataArray <https://docs.xarray.dev/en/stable/generated/xarray.DataArray.html>`_. Output is provided either in a `geopandas.GeoDataFrame <https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.html>`_ or in an `xarray.DataArray <https://docs.xarray.dev/en/stable/generated/xarray.DataArray.html>`_.

Parts of the data setup functions and of the tracking function are based on the `ConTrack - Contour Tracking <https://github.com/steidani/ConTrack>`_ tool developed by `Daniel Steinfeld <https://github.com/steidani>`_. 

**Important information:**

* The package is still under construction and therefore, major errors can occur. 
* Free software: MIT license
* Further documentation about the implemented methods can be found in my `master thesis <https://occrdata.unibe.ch/students/theses/msc/406.pdf>`_

.. end_intro

.. start_installation

Installation
-------------

Stable release
~~~~~~~~~~~~~~~
To install WaveBreaking, run this command in your terminal:
 
..  code-block:: 

        pip install wavebreaking

This is the preferred method to install WaveBreaking, as it will always install the most recent stable release. 
Your virtual environment is automatically checked for the necessary dependencies. 
After the installation, you can start calculating RWB events by following the tutorial below.

From sources
~~~~~~~~~~~~~

The sources for WaveBreaking can be downloaded in two different ways. You can either install WaveBreaking directly from the GitHub repository:

..  code-block:: 

        pip install git+https://github.com/skaderli/WaveBreaking

Or you can clone the GitHub repository first and then install WaveBreaking locally. For that, start with setting the working directory and cloning the repository.

..  code-block:: 

        git clone https://github.com/skaderli/WaveBreaking.git
        cd /path/to/local/WaveBreaking

Second, set up the conda environment and install the necessary dependencies (this may take some time):

..  code-block:: 

        conda create -y -n wb_env
        conda env update -f environment.yml -n wb_env

Now activate the environment and install the WaveBreaking package locally by using the developer mode “-e”:

.. code-block::

        conda activate wb_env
        pip install -e .

To check if the installation was successful, perform some tests:

.. code-block::
 
        python -m unittest tests.test_wavebreaking
        
.. end_installation

.. start_tutorial_part1

Tutorial
---------

This tutorial shows how to calculate RWB events step by step. After successfully installing WaveBreaking, the module needs to be imported. Make sure that the Python kernel with the correct virtual environment (where WaveBreaking is installed) is running.

.. code-block:: python

        import wavebreaking as wb
        
More information about the functions presented below can be found in the `documentation <https://wavebreaking.readthedocs.io/en/latest/modules.html>`_.
   
Data pre-processing:
~~~~~~~~~~~~~~~~~~~~~   

Optionally, the variable intended for the RWB calculations can be smoothed. The smoothing routine applies by default a 5-point smoothing (not diagonally) with a double-weighted center and an adjustable number of smoothing passes. Since the smoothing is based on the scipy.ndimage.convolve function, array-like weights and the mode for handling boundary values can be passed as an argument. This routine returns a xarray.DataArray with the variable "smooth_<variable>". 

.. code-block:: python

        # read data
        import xarray as xr
        demo_data = xr.open_dataset("tests/data/demo_data.nc")

        # smooth variable with 5 passes
        import numpy as np
        smoothed = wb.calculate_smoothed_field(data=demo_data.PV, 
                                               passes=5,
                                               weights=np.array([[0, 1, 0], [1, 2, 1], [0, 1, 0]]), # optional
                                               mode="wrap") # optional
        
The wavebreaking module calculates the intensity for each identified event, if an intensity field is provided. In my master thesis, the intensity is represented by the momentum flux derived from the product of the (daily) zonal deviations of both wind components. The routine creates a xarray.DataArray with the variable "mflux". More information can be found in my `master thesis <https://occrdata.unibe.ch/students/theses/msc/406.pdf>`_.

.. code-block:: python

        # calculate momentum flux
        mflux = wb.calculate_momentum_flux(u=demo_data.U, 
                                           v=demo_data.V)
        
                                   
Contour calculation:
~~~~~~~~~~~~~~~~~~~~
       
All RWB indices are based on a contour line representing the dynamical tropopause. The "calculate_contours()" function calculates the dynamical tropopause on the desired contour levels (commonly the 2 PVU level for Potential Vorticity). The function supports several contour levels at a time which allows for processing data of both hemispheres at the same time (e.g., contour levels -2 and 2). The contour calculation is included in the RWB index functions and doesn't need to be performed beforehand. 

If the input field is periodic, the parameter "periodic_add" can be used to extend the field in the longitudinal direction (default 120 degrees) to correctly extract the contour at the date border. With "original_coordinates = False", array indices are returned (used for the index calculations) instead of original coordinates. The routine returns a geopandas.GeoDataFrame with a geometry column and some properties for each contour. 

.. code-block:: python

        # calculate contours
        contours = wb.calculate_contours(data=smoothed, 
                                         contour_levels=[-2, 2], 
                                         periodic_add=120, # optional
                                         original_coordinates=True) # optional
        

Index calculation:
~~~~~~~~~~~~~~~~~~~

All three RWB indices perform the contour calculation before identifying the RWB events. For the streamer index, the default parameters are taken from `Wernli and Sprenger (2007)`_ (and `Sprenger et al. 2017`_) and for the overturning index from `Barnes and Hartmann (2012)`_. If the intensity is provided (momentum flux, see data pre-processing), it is calculated for each event. All index functions create a geopandas.GeoDataFrame with a geometry column and some properties for each event. 

.. code-block:: python

        # calculate streamers
        streamers = wb.calculate_streamers(data=smoothed, 
                                           contour_levels=[-2, 2], 
                                           geo_dis=800, # optional
                                           cont_dis=1200, # optional
                                           intensity=mflux, # optional
                                           periodic_add=120) # optional
                            
.. code-block:: python                  

        # calculate overturnings
        overturnings = wb.calculate_overturnings(data=smoothed, 
                                                 contour_levels=[-2, 2], 
                                                 range_group=5, # optional
                                                 min_exp=5, # optional
                                                 intensity=mflux, # optional
                                                 periodic_add=120) # optional
        
.. code-block:: python
 
        # calculate cutoffs
        cutoffs = wb.calculate_cutoffs(data=smoothed, 
                                       contour_levels=[-2, 2], 
                                       min_exp=5, # optional
                                       intensity=mflux, # optional
                                       periodic_add=120) # optional
                                       
Event classification:
~~~~~~~~~~~~~~~~~~~~~~

The event classification is based on selecting the events of interest from the geopandas.GeoDataFrame provided by the index calculation functions. 

Some suggested classifications:

.. code-block:: python

        # stratospheric and tropospheric (only for streamers and cutoffs)
        stratospheric = events[events.mean_var >= contour_level]
        tropospheric = events[events.mean_var < contour_level]
        
        # anticyclonic and cyclonic by intensity
        anticyclonic = events[events.intensity >= 0]
        cyclonic = events[events.intensity < 0]
        
        # anticyclonic and cyclonic by orientation (only for overturning events)
        anticyclonic = events[events.orientation == "anticyclonic"]
        cyclonic = events[events.orientation == "cyclonic"]


In addition, a subset of events with certain characteristics can be selected, e.g. the 10% largest events:

.. code-block:: python

        # 10 percent largest events
        large = events[events.area >= events.area.quantile(0.9)]


Transform to DataArray:
~~~~~~~~~~~~~~~~~~~~~~~

To calculate and visualize the occurrence of RWB events, it comes in handy to transform the coordinates of the events into a xarray.DataArray. The "to_xarray" function flags every grid cell where an event is present with the value 1. Before the transformation, it is suggested to classify the events first and only use for example stratospheric events. 

.. code-block:: python

        # classify events
        stratospheric = streamers[streamers.mean_var.abs() >= 2]
        
        # transform to xarray.DataArray
        flag_array = wb.to_xarray(data=smoothed, 
                                  events=stratospheric)

        
Visualization: 
~~~~~~~~~~~~~~~

WaveBreaking provides two options to do a first visual analysis of the output. Both options are based on the xarray.DataArray with the flagged grid cells from the "to_xarray" function. 

To analyze a specific large scale situation, the RWB events on a single time steps can be plotted:

.. code-block:: python

        # import cartopy for projection
        import cartopy.crs as ccrs
        
        wb.plot_step(flag_data=flag_array, 
                     data=smoothed, 
                     step="1959-06-05T06", #index or date
                     contour_level=[-2, 2], # optional
                     proj=ccrs.PlateCarree(), # optional
                     size=(12,8), # optional
                     periodic=True, # optional
                     labels=True,# optional
                     levels=None, # optional
                     cmap="Blues", # optional
                     color_events="gold", # optional
                     title="") # optional

.. end_tutorial_part1

.. image:: https://raw.githubusercontent.com/skaderli/WaveBreaking/main/docs/figures/plot_step.png
    :alt: plot step 
    
.. start_tutorial_part2  
    
The analyze Rossby wave breaking from a climatological perspective, the occurrence (for specific seasons) can be plotted:

.. code-block:: python

        wb.plot_clim(flag_data=flag_array, 
                     seasons=None, # optional
                     proj=ccrs.PlateCarree(), # optional
                     size=(12,8), # optional
                     smooth_passes=0, # optional
                     periodic=True, # optional
                     labels=True, # optional
                     levels=None, # optional
                     cmap=None, # optional
                     title="") # optional

.. end_tutorial_part2

.. image:: https://raw.githubusercontent.com/skaderli/WaveBreaking/main/docs/figures/plot_climatology.png
    :alt: plot climatology 

.. start_tutorial_part3
    
Event tracking:
~~~~~~~~~~~~~~~~

Last but not least, WaveBreaking provides a routine to track events over time. Beside the time range of the temporal tracking, two methods for defining the spatial coherence are available. Events receive the same label if they either spatially overlap (method "by_overlapping") or if the centre of mass lies in a certain radius (method "by_radius"). Again, it is suggested to classify the events first and only use for example stratospheric events. This routine adds a column "label" to the events geopandas.GeoDataFrame.

.. code-block:: python

        # classify events
        anticyclonic = overturnings[overturnings.orientation == "anticyclonic"]

        # track events
        tracked = wb.event_tracking(events=anticyclonic, 
                                    time_range=6, #time range for temporal tracking in hours
                                    method="by_radius", #method for tracking
                                    radius=1000) #radius in km for method "by_radius"

The result can be visualized by plotting the paths of the tracked events:

.. code-block:: python
        
        wb.plot_tracks(data=smoothed,
                       events=tracked,  
                       proj=ccrs.PlateCarree(), # optional
                       size=(12,8), # optional
                       min_path=0, # optional
                       plot_events=True, # optional
                       labels=True, # optional
                       title="") # optional
                       
                  
.. end_tutorial_part3
 
.. image:: https://raw.githubusercontent.com/skaderli/WaveBreaking/main/docs/figures/plot_tracks.png
    :alt: plot tracks

Credits
---------

* The installation guide is to some extend based on the `ConTrack - Contour Tracking <https://github.com/steidani/ConTrack>`_ tool developed by `Daniel Steinfeld <https://github.com/steidani>`_. 

* This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
