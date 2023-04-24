"""
This file is part of WaveBreaking.

WaveBreaking provides indices to detect, classify and track Rossby Wave Breaking (RWB) in climate and weather data.
The tool was developed during my master thesis at the University of Bern.
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---

Spatial pre-processing functions
"""

__author__ = "Severin Kaderli"
__license__ = "MIT"
__email__ = "severin.kaderli@unibe.ch"

# import modules
import xarray as xr
import wrf as wrf

from wavebreaking.utils.data_utils import check_argument_types, get_dimension_attributes


@check_argument_types(["u", "v"], [xr.DataArray, xr.DataArray])
@get_dimension_attributes("u")
def calculate_momentum_flux(u, v, *args, **kwargs):

    """
    Calculate the momentum flux derived from the product of the deviations of both wind components from the zonal mean.

    Parameters
    ----------
        u : xarray.DataArray
            zonal (x) component of the wind
        v : xarray.DataArray
            meridional (y) component of the wind

    Returns
    -------
        momentum flux: xarray.DataArray
            Data containing the momentum flux
    """

    # calculate deviation from the zonal mean
    u_prime = u - u.mean(kwargs["lon_name"])
    v_prime = v - v.mean(kwargs["lon_name"])

    # mflux is given by the product of the deviation of both wind components
    mflux = u_prime * v_prime
    mflux.name = "mflux"

    return mflux


@check_argument_types(["data"], [xr.DataArray])
@get_dimension_attributes("data")
def calculate_smoothed_field(data, passes, *args, **kwargs):

    """
    Calculate smoothed field based on a 5-point smoothing with double-weighted centre and multiple smoothing passes.
    The smoothing routine is based on the wrf.smooth2d function.

    Parameters
    ----------
        data : xarray.DataArray
            data to smooth
        passes : int or float
            number of smoothing passes of the 5-point smoothing

    Returns
    -------
        smoothed data: xarray.DataArray
            Data containing the smoothed field
    """

    # create wrap around the data to get a correct smoothing at the borders
    smoothed = data.pad({kwargs["lon_name"]: 5}, mode="wrap")

    # perform smoothing and remove wrap
    smoothed = wrf.smooth2d(smoothed, passes).assign_attrs(data.attrs).isel({kwargs["lon_name"]: slice(5, -5)})

    return smoothed
