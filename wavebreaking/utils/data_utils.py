"""
This file is part of WaveBreaking.

WaveBreaking provides indices to detect, classify
and track Rossby Wave Breaking (RWB) in climate and weather data.
The tool was developed during my master thesis at the University of Bern
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---

Check input data and arguments
"""

__author__ = "Severin Kaderli"
__license__ = "MIT"
__email__ = "severin.kaderli@unibe.ch"

# import modules
import numpy as np
import geopandas as gpd
import functools

import warnings

warnings.filterwarnings("ignore")


def check_argument_types(arguments, types):
    """
    decorator to check the type of function arguments
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            for (arg_index, arg_name), arg_type in zip(enumerate(arguments), types):
                if arg_name in kwargs:
                    if not isinstance(kwargs[arg_name], arg_type):
                        errmsg = arg_name + " has to be a " + str(arg_type)[8:-2] + "!"
                        raise TypeError(errmsg)
                else:
                    if not isinstance(args[arg_index], arg_type):
                        errmsg = arg_name + " has to be a " + str(arg_type)[8:-2] + "!"
                        raise TypeError(errmsg)

            return func(*args, **kwargs)

        return wrapper

    return decorator


def check_empty_dataframes(func):
    """
    decorator to check if there is an empty DataFrame
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        for item in args:
            if isinstance(item, gpd.GeoDataFrame):
                if item.empty:
                    errmsg = "geopandas.GeoDataFrame is empty!"
                    raise ValueError(errmsg)
        for key, item in kwargs.items():
            if isinstance(item, gpd.GeoDataFrame):
                if item.empty:
                    errmsg = key + " geopandas.GeoDataFrame is empty!"
                    raise ValueError(errmsg)

        return func(*args, **kwargs)

    return wrapper


def get_time_name(data):
    """
    check for 'time' dimension and return name
    """

    for dim in data.dims:
        if (
            ("units" in data[dim].attrs and "since" in data[dim].attrs["units"])
            or (
                "units" in data[dim].encoding and "since" in data[dim].encoding["units"]
            )
            or (data[dim].dtype == np.dtype("datetime64[ns]"))
            or (dim in ["time"])
        ):
            return dim

    errmsg = "'time' dimension (dtype='datetime64[ns]') not found."
    hint = " Add time dimension with xarray.DataArray.expand_dims('time')."
    raise ValueError(errmsg + hint)


def get_lon_name(data):
    """
    check for 'longitude' dimension and return name
    """

    for dim in data.dims:
        if (
            "units" in data[dim].attrs
            and data[dim].attrs["units"] in ["degree_east", "degrees_east"]
        ) or dim in ["lon", "longitude", "x"]:
            return dim

    errmsg = "'longitude' dimension (units='degrees_east') not found."
    raise ValueError(errmsg)


def get_lat_name(data):
    """
    check for 'latitude' dimension and return name
    """

    for dim in data.dims:
        if (
            "units" in data[dim].attrs
            and data[dim].attrs["units"] in ["degree_north", "degrees_north"]
        ) or dim in ["lat", "latitude", "y"]:
            return dim

    errmsg = "latitude' dimension (units='degrees_north') not found."
    raise ValueError(errmsg)


def get_spatial_resolution(data, dim):
    """
    check resolution of the longitude and latitude coordinate
    """

    delta = abs(np.unique((data[dim].data[1:] - data[dim].data[:-1])))

    if len(delta) > 1:
        errmsg = "No regular grid found for dimension {}.".format(dim)
        raise ValueError(errmsg)

    elif delta[0] == 0:
        errmsg = "Two equivalent coordinates found for dimension {}.".format(dim)
        raise ValueError(errmsg)

    return delta[0]


def get_dimension_attributes(arg_name):
    """
    decorator to get the dimension, size and resolution of the input data
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if arg_name in kwargs:
                data = kwargs[arg_name]
            else:
                data = args[0]

            names = ["time_name", "lon_name", "lat_name"]
            sizes = ["ntime", "nlon", "nlat"]
            resolutions = ["dlon", "dlat"]

            get_dims = [get_time_name, get_lon_name, get_lat_name]

            for name, get_dim in zip(names, get_dims):
                if name not in kwargs:
                    kwargs[name] = get_dim(data)

            for size, name in zip(sizes, names):
                if size not in kwargs:
                    kwargs[size] = len(data[kwargs[name]])

            for res, name in zip(resolutions, names[1:]):
                if res not in kwargs:
                    kwargs[res] = get_spatial_resolution(data, kwargs[name])

            if len(data.dims) > 3:
                err_dims = [
                    dim
                    for dim in data.dims
                    if dim not in [kwargs[name] for name in names]
                ]
                errmsg = "Unexpected dimension(s): {}. Select dimensions first.".format(
                    err_dims
                )
                raise ValueError(errmsg)

            return func(*args, **kwargs)

        return wrapper

    return decorator
