"""
This file is part of WaveBreaking. 

WaveBreaking provides indices to detect, classify and track Rossby Wave Breaking (RWB) in climate and weather data. 
The tool was developed during my master thesis at the University of Bern.
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---

Plotting functions
"""

__author__ = "Severin Kaderli"
__license__ = "MIT"
__email__ = "severin.kaderli@unibe.ch"

# import modules
import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from wavebreaking.utils.data_utils import *
from wavebreaking.utils import plot_utils
from wavebreaking.processing import spatial

@check_argument_types(["flag_data"], [xr.DataArray])
@get_dimension_attributes("flag_data")
def plot_clim(flag_data,
              seasons = None, 
              proj = None, 
              size = (12,8), 
              smooth_passes = 5, 
              periodic = True, 
              labels = True, 
              levels = None, 
              cmap = None, 
              title = "",
              *args,
              **kwargs):

    """
    Creates a simple climatological plot showing the occurrence frequency of the detected events. 

    Parameters
    ----------
        flag_data : xarray.DataArray
            data containing the locations of the events flagged with the value 1
        seasons : list or array, optional
            months of the seasons of which the occurrence frequency should be calculated (e.g. [11,12,1])
        proj : cartopy.crs, optional
            cartopy projection
        size : tuple of integers, optional
            size of the figure in the format (width, height)
        smooth_passes : int or float, optional
            number of smoothing passes of the 5-point smoothing of the occurrence frequencies
        periodic : bool, optional
            If True, the first longitude is added at the end to close the gap in a polar projection
        labels : bool, optional
            If False, no labels are added to the plot
        levels : list or array, optional
            Colorbar levles. If not provided, default levels are used.
        cmap : string, optional
            Name of a valid cmap. If not provided, a default cmap is used. 
        title : string, optional
            Title of the plot

    Returns
    -------
        plot : matplotlib.pyplot
            Climatological plot of the occurrence frequencies.     
    """
    
    #define cartopy projections
    data_crs = ccrs.PlateCarree()
    if proj is None:
        proj = data_crs

    #initialize figure
    fig, ax = plt.subplots(1,1, subplot_kw=dict(projection=proj), figsize=size)

    #calculate occurrence frequencies, if provided for seasons
    if seasons is None:
        freq = xr.where(flag_data>0,1,0).sum(dim=kwargs["time_name"])/kwargs["ntime"]*100
    else:
        ds_season = flag_data.sel({self._time_name:flag_data[kwargs["time_name"]].dt.month.isin(seasons)})
        freq = xr.where(ds_season>0,1,0).sum(dim=kwargs["time_name"])/len(ds_season[kwargs["time_name"]])*100
    
    #perform smoothing
    freq = spatial.calculate_smoothed_field(freq.expand_dims("time"), smooth_passes).isel(time=0)

    #add longitude to ensure that there is no gap in a periodic field
    if periodic == True:
        freq = xr.concat([freq, freq.isel({kwargs["lon_name"]:0}).assign_coords({kwargs["lon_name"]:freq[kwargs["lon_name"]].max()+1})], dim = kwargs["lon_name"])
        
    #define levels
    if levels is None:
        levels = plot_utils.get_levels(freq.max())

    #define cmap
    if cmap is None:
        cmap = plot_utils.get_new_cmap("RdYlBu_r")
    
    #plot frequencies
    p = freq.where(freq>0,-999).plot.contourf(ax=ax, cmap = cmap, levels = levels, transform=data_crs, add_colorbar = False, extend = "max")
    
    #define colorbar
    cax = fig.add_axes([ax.get_position().x1+0.05,ax.get_position().y0,0.015,ax.get_position().height])
    cbar = plot_utils.add_colorbar(p, cax, levels, label = "Occurrence frequency in %")
    
    #add coast lines and grid lines
    ax.add_feature(cfeature.COASTLINE, color="dimgrey")
    ax.gridlines(draw_labels=False, color="black", linestyle="dotted", linewidth = 1.1)

    #plot labels
    if labels == True:
        plot_utils.add_grid_lines(ax)

    #make a circular cut out for the NorthPolarStereo projection
    if proj == ccrs.NorthPolarStereo():
        plot_utils.add_circular_boundary(ax)
        plot_utils.add_circular_patch(ax)
        
    #set title
    ax.set_title(title, fontweight='bold',fontsize=20)

    
@check_argument_types(["flag_data", "data"], [xr.DataArray, xr.DataArray])  
@get_dimension_attributes("flag_data")
def plot_step(flag_data, 
              data, 
              step, 
              contour_level, 
              proj = None,
              size = (12,8), 
              periodic = True,
              labels = True, 
              levels = None, 
              cmap = "Blues", 
              color_events = "gold",
              title = "",
              *args,
              **kwargs):
        
    """
    Creates a plot showing the events at one time step.

    Parameters
    ----------
        flag_data : xarray.DataArray
            Data containing the locations of the events flagged with the value 1
        data : xarray.DataArray
            Data that has been used to calculate the contours and the indices
        step : int or string
            index or name of a time step in the xarray.Dataset
        contour_level : list
            list of contour levels that are shown in the plot
        proj : cartopy.crs, optional
            cartopy projection
        size : tuple of integers, optional
            size of the figure in the format (width, height)
        periodic : bool, optional
            If True, the first longitude is added at the end to close the gap in a polar projection
        labels : bool, optional
            If False, no labels are added to the plot
        levels : list or array, optional
            Colorbar levles. If not provided, default levels are used.
        cmap : string, optional
            Name of a valid cmap. If not provided, a default cmap is used. 
        color_events : string, optional
            Color of the plottet events
        title : string, optional
            Title of the plot

    Returns
    -------
        plot : matplotlib.pyplot
            Plot of one time step.      
    """

    #define cartopy projections
    data_crs = ccrs.PlateCarree()
    if proj is None:
        proj = data_crs

    #initialize figure
    fig, ax = plt.subplots(1,1, subplot_kw=dict(projection=proj), figsize=size)

    #select data
    if (type(step) is str or type(step) is np.datetime64):
        ds = data.sel({kwargs["time_name"]:step}).to_dataset()
        ds["flag"] = flag_data.sel({kwargs["time_name"]:step})
        date = pd.to_datetime(str(ds[kwargs["time_name"]].values)).strftime('%Y-%m-%d')  
    else:
        try: 
            ds = data.isel({kwargs["time_name"]:step}).to_dataset()
            ds["flag"] = flag_data.isel({kwargs["time_name"]:step})
            date = pd.to_datetime(str(ds[kwargs["time_name"]].values)).strftime('%Y-%m-%d')
        except ValueError:
            print('{} as input for step is not supported! Step must be either an index or an coordinate of your data!'.format(step))
            
    #get variable name
    variable = data.name
    
    #add longitude to ensure that there is no gap in a periodic field
    if periodic == True:
        ds = xr.concat([ds, ds.isel({kwargs["lon_name"]:0}).assign_coords({kwargs["lon_name"]:ds[kwargs["lon_name"]].max()+1})], dim = kwargs["lon_name"])
            
    #define levels
    if levels is None:
        levels = plot_utils.get_levels(data.max())
            
    #plot variable field, contour line and flag_variable
    p = ds[variable].plot.contourf(ax=ax, cmap = cmap, levels = levels, transform=data_crs, add_colorbar = False, extend = "max", alpha = 0.8)
    c = ds[variable].plot.contour(ax=ax, levels = np.asarray(contour_level).tolist(), transform=data_crs, linewidths = 2, colors = "black")
    f = ds["flag"].where(ds["flag"]>0).plot.contourf(ax=ax, colors = ["white", color_events], levels = [0,0.5], transform=data_crs, add_colorbar = False)
    
    #add the date to the figure
    plt.text(0.99,0.99, "Date: " + str(date), fontsize = 12, fontweight = "bold", ha='right', va='top', transform=ax.transAxes)
    
    #define colorbar
    if all(x in data.attrs for x in ["units", "long_name"]):
        cbar_label = data.long_name + " [" +data.units + "]"
    else:
        cbar_label = None
    cax = fig.add_axes([ax.get_position().x1+0.05,ax.get_position().y0,0.015,ax.get_position().height])
    cbar = plot_utils.add_colorbar(p, cax, levels, label = cbar_label)

    #add coast lines and grid lines
    ax.add_feature(cfeature.COASTLINE, color="dimgrey")
    ax.gridlines(draw_labels=False, color="black", linestyle="dotted", linewidth = 1.1)

    #plot labels
    if labels == True:
        plot_utils.add_grid_lines(ax)

    #make a circular cut out for the NorthPolarStereo projection
    if proj == ccrs.NorthPolarStereo():
        plot_utils.add_circular_boundary(ax)
        plot_utils.add_circular_patch(ax)
        
    #set title
    ax.set_title(title, fontweight='bold',fontsize=20)
    
    
@check_argument_types(["data", "events"], [xr.DataArray, gpd.GeoDataFrame])  
@check_empty_dataframes
@get_dimension_attributes("data")
def plot_tracks(data, 
                events,
                proj = None, 
                size = (12,8),
                min_path = 0, 
                plot_events = False, 
                labels = True, 
                title = "",
                *args,
                **kwargs):

    """
    Creates a plot showing the tracks of the temporally coherent events. 

    Parameters
    ----------
        data : xarray.DataArray
            Data that has been used to calculate the contours and the indices
        events: pd.DataFrame
            DataFrame with the date, coordinates and label of the identified events
        proj : cartopy.crs, optional
            cartopy projection
        size : tuple of integers, optional
            size of the figure in the format (width, height)
        min_path: int, optional
            Minimal number of time steps an event has to be tracked
        plot_events: bool, optional
            If True, the events are also plotted by a shaded area
        labels: bool, optional
            If False, no labels are added to the plot
        title: string, optional
            Title of the plot

    Returns
    -------
        plot : matplotlib.pyplot
            Plot of the tracks
    """

    #define cartopy projections
    data_crs = ccrs.PlateCarree()
    if proj is None:
        proj = data_crs

    #initialize figure
    fig, ax = plt.subplots(1,1, subplot_kw=dict(projection=proj), figsize=size)

    #group event data by label and plot each path
    for name, group in events.groupby("label"):
        if len(group)>min_path:

            lons = np.asarray(group.com.tolist())[:,0]
            lats = np.asarray(group.com.tolist())[:,1]

            #plot start point of each path
            ax.scatter(lons[0],lats[0], s=14, zorder=10,facecolors='none', edgecolor='black', transform=data_crs)
            ax.plot(lons[0], lats[0], ".", color = "red", transform=data_crs, alpha =0.7)

            #plot the coordinates of the events
            if plot_events == True:
                group.plot(ax = ax, color = "black", transform = data_crs, alpha = 0.2, markersize = 3)
                
            max_lon = max(data[kwargs["lon_name"]].values)+1
            min_lon = min(data[kwargs["lon_name"]].values)

            #check if the path needs to be split due to a crossing of the date border
            diffs = abs(np.diff(lons))>kwargs["nlon"]/2
            split = [[i-1,i] for i in np.where(diffs)[0]+1]
            no_split = [[i-1,i] for i in np.where(~diffs)[0]+1]

            #plot paths that do not need to be split
            for item in no_split:
                ev_seg = np.asarray(group[item[0]:item[1]+1].com.tolist())

                ax.plot(ev_seg[:,0], ev_seg[:,1], "-", transform=data_crs, color = "black")

            #plot paths that have to be split
            for item in split:
                ev_seg = np.asarray(group[item[0]:item[1]+1].com.tolist())

                #upstream split
                if np.diff(ev_seg[:,1])<0:
                    lons_plot = [(ev_seg[0,0], max_lon), (min_lon, ev_seg[1,0])]
                    lon_diffs = np.diff(lons_plot)

                    m = np.diff(ev_seg[:,1])[0]/np.sum(lon_diffs)
                    lats_plot = [(ev_seg[0,1],ev_seg[0,1]+lon_diffs[0][0]*m),(ev_seg[1,1]-lon_diffs[1][0]*m, ev_seg[1,1])]

                #downstream split
                else:
                    lons_plot = [(ev_seg[0,0], min_lon), (max_lon, ev_seg[1,0])]
                    lon_diffs = np.diff(lons_plot)

                    m = np.diff(ev_seg[:,1])[0]/np.sum(lon_diffs)
                    lats_plot = [(ev_seg[0,1],ev_seg[0,1]+lon_diffs[0][0]*m),(ev_seg[1,1]-lon_diffs[1][0]*m, ev_seg[1,1])]
                
                #plot splitted segments
                for lon, lat in zip(lons_plot, lats_plot):
                    ax.plot(lon, lat, "-", "-", transform=data_crs, color = "black")        

    #plot invisible data to get a plot of the full grid       
    p = data.isel({kwargs["time_name"]:0}).plot.contourf(ax=ax, transform=data_crs, add_colorbar = False, alpha =0)

    #add coast lines and grid lines
    ax.add_feature(cfeature.COASTLINE, color="dimgrey")
    ax.gridlines(draw_labels=False, color="black", linestyle="dotted", linewidth = 1.1)

    #plot labels
    if labels == True:
        plot_utils.add_grid_lines(ax)

    #make a circular cut out for the NorthPolarStereo projection
    if proj == ccrs.NorthPolarStereo():
        plot_utils.add_circular_boundary(ax)
        plot_utils.add_circular_patch(ax)
        
    #set title
    ax.set_title(title, fontweight='bold',fontsize=20)
    