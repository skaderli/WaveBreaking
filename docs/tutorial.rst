========
Tutorial
========

This tutorial shows how to calculate Rossby wave breaking events step by step. After successfully installing WaveBreaking, the module needs to be imported. Make sure that the Python kernel with the correct virtual environment (where WaveBreaking is installed) is running.

.. code-block:: python

        import wavebreaking as wb
   
        
Data pre-processing:
~~~~~~~~~~       

Optionally, the variable intended for the wave breaking calculations can be smoothed. The smoothing routine applies a 5-point smoothing (not diagonally) with a double-weighted center and an adjustable number of smoothing passes. This routine returnes a xarray.DataArray with the variable "smooth_variable". 

.. code-block:: python

        #read your data
        import xarray as xr
        demo_data = xr.open_dataset("tests/tests_data/test_data.nc")

        #smooth variable with 5 passes
        smoothed = wb.calculate_smoothed_field(data = demo_data.variable, 
                                                      passes = 5)
        
        
The wavebreaking module can calculate the intensity for each identified event. For that, the intensity field needs to be provided or calculated before the event identification. Here, the momentum flux derived from the product of the (daily) zonal deviations of both wind components is used as the intensity. This routine creates a xarray.DataArray with the variable "mflux". More information can be found in my `master thesis <https://occrdata.unibe.ch/students/theses/msc/406.pdf>`_.

.. code-block:: python

        #calculate momentum flux (wind data not included in the demo data)
        mflux = wb.calculate_momentum_flux(u = demo_data.variable, 
                                           v = demo_data.variable)
        
                                   
Contour calculation:
~~~~~~~~~~
       
Both Rossby wave breaking indices are based on a contour line representing the dynamical tropopause. The "calculate_contours()" function calculates the dynamical tropopause on the desired contour levels (commonly the 2 PVU level for Potential Vorticity). The function supports several contour levels at a time which allows for the contour calculation on both Hemispheres at the same time (levels -2 and 2). 

If the input field is periodic, the parameter "periodic_add" can be used to extend the field in the longitudinal direction (default 120 degrees) to correctly extract the contour at the date border. With "original_coordinates = False", array indices are returned (used for the index calculations), instead of original coordinates. This routines returns a geopandas.GeoDataFrame with a geometry column and some properties for each contour. 

.. code-block:: python

        #calculate contours
        contours = wb.calculate_contours(data = smoothed, 
                                         contour_levels = [-2, 2], 
                                         periodic_add = 120, 
                                         original_coordinates = True)
        

Index calculation:
~~~~~~~~~~

All three RWB indices perform the contour calculation before identifying the RWB events. For the streamer index, the default parameters are taken from `Wernli and Sprenger (2007)`_ (and `Sprenger et al. 2017`_) and for the overturning index from `Barnes and Hartmann (2012)`_. If the intensity is provided (momentum flux, see data pre-processing), it is calculated for each event. All index functions create a geopandas.GeoDataFrame with a geometry column and some properties for each event. 

.. code-block:: python

        #calculate streamers
        streamers = wb.calculate_streamers(data = smoothed, 
                                           contour_levels = [-2, 2], 
                                           geo_dis = 800,
                                           cont_dis = 1200,
                                           intensity = mflux,
                                           periodic_add = 120)
                            
.. code-block:: python                  

        #calculate overturnings
        overturnings = wb.calculate_overturnings(data = smoothed, 
                                                 contour_levels = [-2, 2], 
                                                 range_group = 500, 
                                                 min_exp = 5, 
                                                 intensity = mflux,
                                                 periodic_add = 120)
        
.. code-block:: python
 
        #calculate cutoffs
        cutoffs = wb.calculate_cutoffs(data = smoothed, 
                                       contour_levels = [-2, 2], 
                                       min_exp = 5,
                                       intensity = mflux, 
                                       periodic_add = 120)


Transform to xarray.DataArray:
~~~~~~~~~~

To calculate and visualize the occurrence of Rossby wave breaking, it comes in handy to transform the coordinates of the events into a xarray.DataArray. The "to_xarray" function flags every grid cell where an event is present with the value 1. Before the transformation, it is suggested to filter the geopandas.GeoDataFrame for the desired events (e.g., stratospheric events with Potential Vorticity values larger than 2 PVU).

.. code-block:: python

        #filter events
        f_events = streamers[streamers.mean_var >= 2]
        
        #transform to xarray.DataArray
        flag_array = wb.to_xarray(data = smoothed, #data used for the index calculation (to receive the same dimensions)
                                  events = f_events)

        
Visualization: 
~~~~~~~~~~

WaveBreaking provides two options to do a first visual analysis of the output. Both options are based on the xarray.DataArray with the flagged grid cells from the "to_xarray" function. 

To analyze a specific large scale situation, the wave breaking events on a single time steps can be plotted:

.. code-block:: python

        #import cartopy for projection
        import cartopy.crs as ccrs
        
        wb.plot_step(flag_data = flag_array, 
                     data = smoothed, 
                     step = "1979-06-18", #index or date (this date is not in the demo data)
                     contour_level = [2],
                     proj = ccrs.NorthPolarStereo(), #cartopy projection, optional
                     size = (12,8), 
                     periodic = True, 
                     labels = True,
                     levels = None, 
                     cmap = "Blues",
                     color_events = "gold", 
                     title = "")

.. image:: docs/figures/plot_step.png
    :alt: plot step 
    
The analyze Rossby wave breaking from a climatological perspective, the occurrence (for specific seasons) can be plotted:

.. code-block:: python

        wb.plot_clim(flag_data = flag_array, 
                     seasons = None,
                     proj = ccrs.NorthPolarStereo(), #cartopy projection, optional
                     size = (12,8), 
                     smooth_passes = 0,
                     periodic = True, 
                     labels = True,
                     levels = None, 
                     cmap = None, 
                     title = "")

.. image:: docs/figures/plot_climatology.png
    :alt: plot climatology 
    
Event tracking:
~~~~~~~~~~~

Last but not least, the wave breaking module provides a routine to track events over time. Events that overlap between two time steps receive the same label. Again, it is suggested to filter the events first. This routine adds a column "label" to the events geopandas.GeoDataFrame.

.. code-block:: python

        #filter events
        f_events = streamers[streamers.mean_var >= 2][::2] #use every second event for clarity

        #track events
        wb.event_tracking(events = f_events, 
                          time_range = 24) #time range for temporal tracking in hours

The result can be visualized by plotting the paths of the tracked events:

.. code-block:: python
        
        wb.plot_tracks(data = smoothed, #data used for the index calculation 
                       events = f_events,  
                       proj = ccrs.NorthPolarStereo(), #cartopy projection, optional
                       size = (12,8),
                       min_path = 0, #minimal number of paths
                       plot_events = False, #plot events as grey shaded area
                       labels = True,
                       title = "")
                       
.. image:: docs/figures/plot_tracks.png
    :alt: plot tracks

