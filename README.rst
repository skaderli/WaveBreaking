.. image:: https://img.shields.io/pypi/v/wavebreaking.svg
        :target: https://pypi.python.org/pypi/wavebreaking

============
WaveBreaking - Detection, Classification and Tracking of Rossby Wave Breaking
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

Simple instructions
~~~~~~~~
By using pip, WaveBreaking can be directly installed in a virtual environment:
 
..  code-block:: 

        pip install wavebreaking ### NOT AVAILABLE YET

The evironment is automatically checked for the necessary dependencies. After the installation you can start calculating wave breaking events by following the tutorial below. 

The WaveBreaking GitHub repository can also be directly installed by using pip:

..  code-block:: 

        pip install git+https://github.com/skaderli/WaveBreaking

Advanced instructions
~~~~~~~~
For developers, it is suggested to clone the WaveBreaking repository first and then install the tool locally. After setting the working directory, cloning the repository creates a directory named WaveBreaking. 

..  code-block:: 

        cd /path/to/local/workspace
        git clone https://github.com/skaderli/WaveBreaking.git

In the next step, the conda environment is set up and the necessary dependencies are installed:

..  code-block:: 

        conda create -y -n wb_dev
        conda env update -f environment_dev.yml -n wb_dev

Now the environment can be activated and the WaveBreaking package can be locally installed by using the developer mode "-e":

.. code-block::

        conda activate wb_dev
        pip install -e .

To check if the installation was sucessful, some tests can be performed:

.. code-block::
 
        python -m pytest
        

Tutorial
--------

This tutorial shows how to calculate Rossby wave breaking events step by step. After successfully installling the wavebreaking package, the wavebreaking module needs to be imported. Make sure that the Python kernel with the correct virtual environment (where the wavebreking package is instaled) is running.

.. code-block:: python

        from wavebreaking import wavebreaking
        
Read data:
~~~~~~~~~~

Input data is only accepted in a NetCDF-file with two spatial and one temporal dimensions. There are two options to read data: Either directly as a NetCDF-file or as a xarray.DataSet: 

.. code-block:: python

        #input data 
        import xarray as xr
        file = "tests/data/test_data.nc"
        ds = xr.open_dataset(file)

        #initiate wavebreaking class and read data
        wb = wavebreaking(file) #or
        wb = wavebreaking(ds)
        
        #data can also be read in explicitly
        wb = wavebreaking()
        wb.read(file) #or
        wb.read_xarray(ds)
        
Data pre-processing:
~~~~~~~~~~       

Optionally, the variable intended for the wave breaking calculations can be smoothed. The smoothing routine applies a 5-point smoothing (not diagonally) with a double-weighted center and an adjustable number of smoothing passes. This routine creates a xarray.DataArray with the variable "smooth_variable". 

.. code-block:: python

        #smooth variable with 5 passes
        wb.calculate_smoothed_field("variable", passes = 5)
        
        #access xarray.DataArray
        wb["smooth_variable"]
        
The wavebreaking module can calculate the intensity for each identified breaking event. For that, the intensity field needs to be calculated before the event identification. Here, the momentum flux is calculated as the product of the (daily) zonal deviation of both wind components. This routine creates a xarray.DataArray with the variable "mflux". More information can be found in my `master thesis <https://occrdata.unibe.ch/students/theses/msc/406.pdf>`_. If the momentum flux is not calculated, the intensity of the events is not provided.

.. code-block:: python

        #calculate momentum flux
        wb.calculate_momentum_flux(variable_zonal = "zonal", variable_meridional = "meridional", dtime = "1D")
        
        #access xarray.DataArray
        wb["mflux"]
                                   
Contour calculation:
~~~~~~~~~~
       
Both Rossby wave breaking indices are based on a contour line representing the dynamical tropopause. The "get_contours()" function calculates the dynamical tropopause on a specific level (commonly the 2 PVU level for Potential Vorticity). If the input field is periodic, the parameter "periodic_add" can be used to extend the field in the longitudinal direction (default 120 degrees) to correctly extract the contour at the date border. This routines creates a pandas.DataFrame with the coordinates and some properties of the contour line.

.. code-block:: python

        #calculate contours
        wb.get_contours(variable = "smooth_variable", level = 2, periodic_add = 120)
        
        #access pandas.DataFrame
        wb.contours 
        

Index calculation:
~~~~~~~~~~

Now the index calculation can be performed based on the identified contour lines. For the streamer index, the default parameters are taken from `Wernli and Sprenger (2007)`_ (and `Sprenger et al. 2017`_) and for the overturning index from `Barnes and Hartmann (2012)`_. Both index functions create a pd.DataFrame with the coordinates and some properties of the events.

.. code-block:: python

        #calculate events
        wb.get_streamers(geo_dis = 800, cont_dis = 1500)
        wb.get_overturnings(range_group = 500, min_exp = 5)
        
        #access pandas.DataFrame
        wb.streamers
        wb.overturnings

Transform to DataSet:
~~~~~~~~~~

To calculate and visualize the occurrence of Rossby wave breaking, it comes in handy to transform the coordinates of the events into a xarray.DataSet. The "to_xarray" function flags every grid cell where an event is present with the value 1. Before the transformation, it is suggested to filter the pandas.DataFrame for the desired events (e.g., stratospheric events with a Potential Vorticity value larger than 2 PVU).

.. code-block:: python

        #filter events
        f_events = wb.streamers[wb.streamers.mean_var >= 2]
        
        #transform to xarray
        wb.to_xarray(f_events, name = "flag")
        
        #access xarray.DataSet
        wb.flag
        
Visualization: 
~~~~~~~~~~

The wavebreaking module provides two options to do a first visual analysis of the ouput. Both options are based on the xarray.Dataset with the flagged grid cells from the "to_xarray" function. 

To analyze a specific large scale situation, the wave breaking events on a single time steps can be plotted:

.. code-block:: python
        
        wb.plot_step(variable = "smooth_variable", #variable used for contour calculation
                     flag_variable = "stratos_streamers", 
                     contour_level = [2], 
                     step = "1900-01-01", #date or index
                     proj = "NorthPolarStereo", #name of cartopy projection,
                     labels = True, 
                     levels = None, #levels color bar
                     cmap = None, 
                     title = "Potential vorticity streamers"
                     )

.. image:: docs/plot_step.png
    :alt: plot step 
    
The analyze Rossby wave breaking from a climatological perspective, the occurrence (for specific seasons) can be plotted:

.. code-block:: python

        wb.plot_clim(variable = "stratos_streamers", 
                     seasons = None, #specify months  
                     smooth_passes = 2, #of the occurrence frequency field  
                     proj = "NorthPolarStereo", 
                     periodic = True, 
                     labels = True, 
                     levels = None, 
                     cmap = None, 
                     title = "Climatological analysis"
                     )

.. image:: docs/plot_climatology.png
    :alt: plot climatology 
    
Event tracking:
~~~~~~~~~~~

Last but not least, the wave breaking module provides a routine to track events over time. Events, that overlap between two time steps receive the same label. Again, it is suggested to filter the events first. This routine adds a label column to the events pandas.DataFrame.

.. code-block:: python

        #filter events
        f_events = wb.streamers[wb.streamers.mean_var >= 2]

        #track events
        wb.event_tracking(f_events, 
                          box = False #if True, a rectangular box is used for the tracking
                          )

The result can be viszalized by plotting the paths of the tracked events:

.. code-block:: python
        
        wb.plot_tracks(events = wb.labeled_events, 
                       proj = "NorthPolarStereo", 
                       min_path = 0, #minimal number of steps per tracked event
                       labels = True,
                       plot_events = False, #if True, the grid cells of the events are shaded
                       title = "Events paths"
                       )
                       
.. image:: docs/plot_tracks.png
    :alt: plot tracks


Credits
-------

* The installation giude above is to some extend based on the `ConTrack - Contour Tracking <https://github.com/steidani/ConTrack>`_ tool developed by `Daniel Steinfeld <https://github.com/steidani>`_. 

* This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
