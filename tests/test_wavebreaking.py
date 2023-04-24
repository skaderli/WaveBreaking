"""
This file is part of WaveBreaking. 

WaveBreaking provides indices to detect, classify and track Rossby Wave Breaking (RWB) in climate and weather data. 
The tool was developed during my master thesis at the University of Bern
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---
Test for WaveBreaking
"""

__author__ = "Severin Kaderli"
__license__ = "MIT"
__email__ = "severin.kaderli@unibe.ch"

# import modules
import unittest
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
from shapely.geometry import MultiPoint
import functools
import logging

# import and silent tqdm progress bar
import tqdm
def silent_tqdm(it, *args, **kwargs):
    return it
tqdm.tqdm = silent_tqdm

from wavebreaking.utils.data_utils import *
from wavebreaking.utils.index_utils import *
from wavebreaking.utils.plot_utils import *

from wavebreaking.processing.spatial import *
from wavebreaking.processing.events import *

from wavebreaking.indices.contour_index import *
from wavebreaking.indices.streamer_index import *
from wavebreaking.indices.overturning_index import *
from wavebreaking.indices.cutoff_index import *

data = xr.open_dataset('tests/tests_data/test_data.nc')

# create decorater to disable logging for testing
def disablelogging(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        logger = logging.getLogger('wavebreaking')
        previousloglevel = logger.getEffectiveLevel()
        logger.setLevel(level=logging.WARNING)
        try:
            return func(*args, **kwargs)
        finally:
            logger.setLevel(previousloglevel)
    return wrapper


class test_data_utils(unittest.TestCase):
    
    def test_check_argument_types(self):
        
        @check_argument_types(["data"], [xr.DataArray])
        def to_be_decorated(data, *args, **kwargs):
            return None 
            
        self.assertRaises(TypeError, lambda: to_be_decorated(""))
        
    def test_get_dimension_attributes(self):
        
        @get_dimension_attributes("data")
        def to_be_decorated(data, *args, **kwargs):
            
            self.assertEqual(kwargs["time_name"], "time")
            self.assertEqual(kwargs["lon_name"], "lon")
            self.assertEqual(kwargs["lat_name"], "lat")
            
            self.assertEqual(kwargs["ntime"], 15)
            self.assertEqual(kwargs["nlon"], 360)
            self.assertEqual(kwargs["nlat"], 90)
            
            self.assertEqual(kwargs["dlon"], 1)
            self.assertEqual(kwargs["dlat"], 1)
            
        to_be_decorated(data.variable)
        
        
class test_index_utils(unittest.TestCase):
    
    def test_add_logger(self):
        
        with self.assertLogs() as captured:

            @add_logger("Logger message")
            def to_be_decorated(*args, **kwargs):
                return None
        
            to_be_decorated()

        self.assertEqual(len(captured.records), 1) # check that there is only one log message
        self.assertEqual(captured.records[0].getMessage(), "Logger message") # and it is the proper one
        
    def test_iterate_time_dimension(self):
        
        @get_dimension_attributes("data")
        @iterate_time_dimension
        def to_be_decorated(data, *args, **kwargs):
            return pd.DataFrame([kwargs["step"].values], columns = ["time"])
        
        df_check = data.time.to_dataframe().reset_index(drop=True)
        self.assertEqual(True, to_be_decorated(data.variable).equals(df_check))
        
    def test_combine_shared(self):
        
        list_in = [[1,2,3], [2,3,4], [5,6]]
        list_out = [[1,2,3,4], [5,6]]
        self.assertEqual(combine_shared(list_in), list_out)
        
    def test_properties_per_timestep(self):
        
        df_in = pd.DataFrame([[10,50]], columns = ["y", "x"])

        df_check = properties_per_timestep("1900-01-01", [df_in], data.variable, data.variable, periodic_add = 120, 
                                           time_name = "time", lat_name = "lat", nlat = 90, nlon = 360, lon_name = "lon")
        
        self.assertIs(type(df_check), gpd.GeoDataFrame)
        self.assertEqual(len(df_check.columns), 6)
       
          
class test_spatial(unittest.TestCase):
    
    def test_calculate_momentum_flux(self):
        
        self.assertIs(type(calculate_momentum_flux(data.variable, data.variable)), xr.DataArray)
        
    def test_calculate_smoothed_filed(self):
        
        self.assertIs(type(calculate_smoothed_field(data.variable, 10)), xr.DataArray)
   
        
class test_events(unittest.TestCase):
    
    def test_to_xarray(self):
        
        date = "1900-01-01"
        name = "test_flag"
        events = gpd.GeoDataFrame(pd.DataFrame([{"date":date}]), geometry = [MultiPoint([(-148.0, 27.0)])])
        
        flag_data = to_xarray(data = data.variable, events = events, name = name)
        
        self.assertIs(type(flag_data), xr.DataArray)
        self.assertEqual(flag_data.sel(lon = -148, lat = 27, time = date).values, 1)
        self.assertEqual(flag_data.name, name)
        
    def test_track_events(self):
        
        date = "1900-01-01"
        events = gpd.GeoDataFrame(pd.DataFrame([{"date":date}]), geometry = [MultiPoint([(-148.0, 27.0)])])
        
        tracked = track_events(events = events, time_range = 24)
        
        self.assertIs(type(tracked), gpd.GeoDataFrame)
        self.assertEqual(tracked.iloc[0].label, 0)
        
        
class test_indices(unittest.TestCase):
    
    @disablelogging
    def test_contour_index(self):
        
        contours_coords = calculate_contours(data = data.variable, 
                                                           contour_level = 2, 
                                                           periodic_add = 120, 
                                                           original_coordinates = True)
        
        contours_index = calculate_contours(data = data.variable, 
                                            contour_level = 2, 
                                            periodic_add = 120, 
                                            original_coordinates = False)
        
        self.assertIs(type(contours_coords), gpd.GeoDataFrame)
        self.assertIs(type(contours_index), gpd.GeoDataFrame)
        
        self.assertEqual(min([p.x for p in contours_coords.iloc[0].geometry.geoms]), -180)
        self.assertEqual(min([p.x for p in contours_index.iloc[0].geometry.geoms]), 0)
        
        cols = ['date', 'exp_lon', 'mean_lat', 'geometry']
        self.assertEqual(contours_coords.columns.to_list(), cols)
        self.assertEqual(contours_index.columns.to_list(), cols)
    
    @disablelogging
    def test_decorator_contour_calculation(self):
        
        @decorator_contour_calculation
        def to_be_decorated(*args, **kwargs):
            contours = kwargs["contours"]
            
            self.assertIs(type(contours), gpd.GeoDataFrame)
            self.assertEqual(min([p.x for p in contours.iloc[0].geometry.geoms]), 0)
            self.assertEqual(contours.columns.to_list(), ['date', 'exp_lon', 'mean_lat', 'geometry'])
            
        to_be_decorated(data.variable, contour_level = 2)
    
    @disablelogging
    def test_streamer_index(self):
        
        streamers = calculate_streamers(data = data.variable, 
                                        contour_level = 2, 
                                        geo_dis = 800,
                                        cont_dis = 1200,
                                        intensity = data.variable,
                                        periodic_add = 120)
        
        self.assertIs(type(streamers), gpd.GeoDataFrame)
        self.assertEqual(len(streamers), 20)
        self.assertEqual(streamers.columns.to_list(), ['date', 'com', 'mean_var', 'area', 'intensity', 'geometry'])
        
    @disablelogging
    def test_overturning_index(self):
        
        overturnings = calculate_overturnings(data = data.variable, 
                                             contour_level = 2, 
                                             range_group = 500, 
                                             min_exp = 5, 
                                             intensity = data.variable,
                                             periodic_add = 120)
        
        self.assertIs(type(overturnings), gpd.GeoDataFrame)
        self.assertEqual(len(overturnings), 12)
        self.assertEqual(overturnings.columns.to_list(), ['date', 'com', 'mean_var', 'area', 'intensity', 'geometry'])
        
    @disablelogging
    def test_cutoff_index(self):
        
        cutoffs = calculate_cutoffs(data = data.variable, 
                                    contour_level = 2,
                                    min_exp = 5,
                                    intensity = data.variable, 
                                    periodic_add = 120)
        
        self.assertIs(type(cutoffs), gpd.GeoDataFrame)
        self.assertEqual(len(cutoffs), 10)
        self.assertEqual(cutoffs.columns.to_list(), ['date', 'com', 'mean_var', 'area', 'intensity', 'geometry'])

# Execute Tests      
if __name__ == '__main__':
    unittest.main()
        