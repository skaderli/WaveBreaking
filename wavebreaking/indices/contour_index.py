"""
This file is part of WaveBreaking.

WaveBreaking provides indices to detect, classify and track Rossby Wave Breaking (RWB) in climate and weather data.
The tool was developed during my master thesis at the University of Bern.
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---

Contour calculation
"""

__author__ = "Severin Kaderli"
__license__ = "MIT"
__email__ = "severin.kaderli@unibe.ch"

# import modules
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import itertools as itertools
from shapely.geometry import MultiPoint
from skimage import measure
import functools

from wavebreaking.utils.data_utils import (
    get_dimension_attributes,
    add_logger)
from wavebreaking.utils.index_utils import (
    iterate_time_dimension)


@get_dimension_attributes("data")
@add_logger("Calculating contours...")
@iterate_time_dimension
def calculate_contours(data,
                       contour_level,
                       periodic_add=120,
                       original_coordinates=True,
                       *args,
                       **kwargs):

    """
    Calculate contour lines of a specific level at each time step. The calculations are based on the measure.find_contours module.
    If periodic_add is provided, the data array is expanded in the longitudinal direction and undulations at the date border are correctly captured.

    Parameters
    ----------
        data : xarray.DataArray
            data for the contour calculation
        contour_level : int or float
            level of the contours
        periodic_add: int or float, optional
            number of longitudes in degree to expand the dataset to correctly capture undulations at the date border
            if the input field is not periodic, use periodic_add = 0
        original_coordinates: bool, optional
            if False, the array indices of the contour lines are returned instead of the original coordinates

    Returns
    -------
        contours: geopandas.GeoDataFrame
            GeoDataFrame containing different characteristics of the contours
                "date": date of the contour line
                "exp_lon": expansion of the contours in the longitudinal direction
                "mean_lat": mean latitude of the contours
                "geometry": MultiPoint object with the coordinate points of the contours in the format (x,y)
    """

    # select variable and time step for the contour calculation
    ds = data.sel({kwargs["time_name"]: kwargs["step"]})
    date = pd.to_datetime(ds.time.values).strftime('%Y-%m-%dT%H')

    # expand field for periodicity
    ds = xr.concat([ds, ds.isel({kwargs["lon_name"]: slice(0, periodic_add)})], dim=kwargs["lon_name"])

    # get contours (indices in array coordinates)
    contours_from_measure = measure.find_contours(ds.values, contour_level)

    # contours in indices for wave breaking calculation
    contours_index_expanded = [np.asarray(list(dict.fromkeys(map(tuple, np.round(item).astype("int"))))) for item in contours_from_measure]

    def contour_to_dataframe(list_of_arrays):

        """
        Calculate different characteristics of the contour line.
        """

        exp_lon = [len(set(item[:, 1])) for item in list_of_arrays]
        mean_lat = [np.round(item[:, 0].mean(), 2) for item in list_of_arrays]
        geo_mp = [MultiPoint(coords[:, ::-1]) for coords in list_of_arrays]

        gdf = gpd.GeoDataFrame(pd.DataFrame({"date": date, "exp_lon": exp_lon, "mean_lat": mean_lat}), geometry=geo_mp)

        return gdf

    if original_coordinates is False:
        # return contours in index coordinates as a geopandas.GeoDataFrame
        return contour_to_dataframe(contours_index_expanded)

    else:

        def map_indices_to_original_grid(list_of_arrays):

            """
            Map contour indices to original input grid.
            """

            # get original coordinates
            contours_index_original = [np.c_[item[:, 0] % kwargs["nlat"], item[:, 1] % kwargs["nlon"]] for item in list_of_arrays]
            contours_index_original = [list(dict.fromkeys(map(tuple, item))) for item in contours_index_original]

            # drop duplicates that have been identified due to the periodicity
            index = np.arange(0, len(contours_index_original))
            index_combination = list(itertools.combinations(index, r=2))

            drop = []
            for combination in index_combination:
                if set(contours_index_original[combination[0]]) < set(contours_index_original[combination[1]]):
                    drop.append(combination[0])
                if set(contours_index_original[combination[1]]) < set(contours_index_original[combination[0]]):
                    drop.append(combination[1])
                if set(contours_index_original[combination[1]]) == set(contours_index_original[combination[0]]):
                    drop.append(combination[1])

            contours_index_original = [np.asarray(contours_index_original[i]) for i in index if i not in drop]

            # select the original coordinates from the indices
            return [np.c_[data[kwargs["lat_name"]].values[item[:, 0].astype("int")],
                          data[kwargs["lon_name"]].values[item[:, 1].astype("int")]] for item in contours_index_original]

        # return contours in original cooridnates as a geopandas.GeoDataFrame
        return contour_to_dataframe(map_indices_to_original_grid(contours_index_expanded))


def decorator_contour_calculation(func):

    """
    decorator to wrap the contour calulcation around the index functions
    """

    @functools.wraps(func)
    def wrapper(data, contour_level, periodic_add=120, original_coordinates=False, *args, **kwargs):

        # pass contours to the wrapped function as a key word argument
        kwargs["contours"] = calculate_contours(data, contour_level, periodic_add, original_coordinates)

        return func(data, contour_level, *args, **kwargs)
    return wrapper
