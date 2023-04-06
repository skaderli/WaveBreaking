#!/usr/bin/env python

"""
Tests for `wavebreaking` package.
    
@author: severinkaderli (Severin Kaderli; severin.kaderli@unibe.ch)
Remark : The test functions are based on the ConTrack - Contour Tracking tool developed by Daniel Steinfeld

"""

import pytest

from wavebreaking import wavebreaking

dataset = 'tests/data/test_data.nc'

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
    except IOError as err:
        assert err.args[0] == "Unkown fileformat. Known formats " \
                               "are netcdf."
def test_read_xarray():
    data = xr.open_dataset(dataset)
    wavebreakings = wavebreaking()
    wavebreakings.read_xarray(data)
    assert type(wavebreakings.ds) is xr.Dataset
    
def test_len(wavebreakings):
    assert len(wavebreakings) == 1

def test_ntime(wavebreakings):
    assert wavebreakings.ntime == 11
    
def test_dimensions(wavebreakings):
    assert wavebreakings.dimensions == ['lat', 'lon', 'time']
    
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
