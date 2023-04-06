#!/usr/bin/env python

"""
Tests for `wavebreaking` package.
    
@author: severinkaderli (Severin Kaderli; severin.kaderli@unibe.ch)
Remark : The test functions are based on the ConTrack - Contour Tracking tool developed by Daniel Steinfeld

"""

import numpy as np
import xarray as xr
import pandas as pd

import matplotlib.pyplot as plt

import warnings
import pytest

from wavebreaking import wavebreaking

dataset = 'tests/data/test_data.nc'
wrong_dataset = 'tests/data/test_data.txt'

@pytest.fixture
def wavebreakings():
    wavebreakings = wavebreaking()
    wavebreakings.read(dataset)
    return wavebreakings

def test_init_empty():
    wavebreakings = wavebreaking()
    assert wavebreakings.ds is None

def test_init_netcdf(wavebreakings):
    assert type(wavebreakings.ds) is xr.Dataset

def test_read_netcdf():
    wavebreakings = wavebreaking(dataset)
    assert type(wavebreakings) == wavebreaking
    
def test_read_wrong():
    try:
        wavebreaking(wrong_dataset)
    except ValueError as err:
        assert err.args[0] == "Unkown fileformat. Known formats are netcdf."
        
def test_read_xarray():
    data = xr.open_dataset(dataset)
    wavebreakings = wavebreaking()
    wavebreakings.read_xarray(data)
    assert type(wavebreakings.ds) is xr.Dataset
    
def test_len(wavebreakings):
    assert len(wavebreakings) == 1

def test_ntime(wavebreakings):
    assert wavebreakings.ntime == 15
    
def test_dimensions(wavebreakings):
    assert wavebreakings.dimensions == ['lon', 'lat', 'time']
    
def test_variables(wavebreakings):
    assert wavebreakings.variables == ['variable']
    
def test_set_up_manually(wavebreakings):
    wavebreakings.set_up(time_name='time',
                     longitude_name='lon',
                     latitude_name='lat')    
    assert wavebreakings._time_name == 'time'
    assert wavebreakings._longitude_name == 'lon'
    assert wavebreakings._latitude_name == 'lat'

def test_set_up_automatic(wavebreakings):
    wavebreakings.set_up()   
    assert wavebreakings._time_name == 'time'
    assert wavebreakings._longitude_name == 'lon'
    assert wavebreakings._latitude_name == 'lat'
    
def test_calculate_momentum_flux(wavebreakings):
    wavebreakings.set_up()  
    wavebreakings.calculate_momentum_flux(variable_zonal = 'variable', variable_meridional = 'variable')
    assert "mflux" in wavebreakings.dataset.keys()
    assert type(wavebreakings.dataset.mflux) is xr.DataArray

def test_calculate_smoothed_field(wavebreakings):
    wavebreakings.set_up() 
    wavebreakings.calculate_smoothed_field(variable = 'variable', passes = 10)
    assert "smooth_variable" in wavebreakings.dataset.keys()
    assert type(wavebreakings.dataset.smooth_variable) is xr.DataArray
    
def test_get_contours(wavebreakings):
    wavebreakings.set_up()
    wavebreakings.get_contours(variable = "variable", level = 2)
    assert type(wavebreakings.contours) is pd.DataFrame
    assert list(wavebreakings.contours.columns) == ['date', 'coordinates', 'level', 'exp_lon', 'mean_lat']
    assert type(wavebreakings._contours_wb_calc) is pd.DataFrame
    
def test_get_streamers(wavebreakings):
    wavebreakings.set_up()
    wavebreakings.get_contours(variable = "variable", level = 2)
    wavebreakings.get_streamers()
    assert type(wavebreakings.streamers) is pd.DataFrame
    assert len(wavebreakings.streamers) == 15
    assert list(wavebreakings.streamers.columns) == ['date', 'coordinates', 'mean_var', 'area', 'com']
    
def test_get_overturnings(wavebreakings):
    wavebreakings.set_up()
    wavebreakings.get_contours(variable = "variable", level = 2)
    wavebreakings.get_overturnings()
    assert type(wavebreakings.overturnings) is pd.DataFrame
    assert len(wavebreakings.overturnings) == 12
    assert wavebreakings.overturnings.columns.tolist() == ['date', 'coordinates', 'mean_var', 'area', 'com']
    
def test_to_xarray(wavebreakings):
    wavebreakings.set_up()
    wavebreakings.get_contours(variable = "variable", level = 2)
    wavebreakings.get_overturnings()
    wavebreakings.to_xarray(events = wavebreakings.overturnings)
    assert "flag" in wavebreakings.dataset.keys()
    assert type(wavebreakings.dataset.flag) is xr.DataArray
    
def test_event_tracking(wavebreakings):
    wavebreakings.set_up()
    wavebreakings.get_contours(variable = "variable", level = 2)
    wavebreakings.get_overturnings()
    wavebreakings.event_tracking(events = wavebreakings.overturnings)
    assert type(wavebreakings.labeled_events) is pd.DataFrame
    assert len(set(wavebreakings.labeled_events.label)) == 7
    assert wavebreakings.labeled_events.columns.tolist() == ['date', 'coordinates', 'mean_var', 'area', 'com', 'label']
    
def test_plot_clim(wavebreakings):
    try:
        wavebreakings.set_up()
        wavebreakings.plot_clim(variable = "variable")
    except:
        assert False
        
def test_plot_step(wavebreakings):
    try:
        wavebreakings.set_up()
        wavebreakings.plot_step(flag_variable = "variable", variable = "variable", contour_level = [2], step = 0)
    except:
        assert False    

def test_plot_tracks(wavebreakings):
    try:
        wavebreakings.set_up()
        wavebreakings.get_contours(variable = "variable", level = 2)
        wavebreakings.get_overturnings()
        wavebreakings.event_tracking(events = wavebreakings.overturnings)
        wavebreakings.plot_tracks(events = wavebreakings.labeled_events)
    except:
        assert False    
