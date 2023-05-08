""""""
"""
This file is part of WaveBreaking.

WaveBreaking provides indices to detect, classify
and track Rossby Wave Breaking (RWB) in climate and weather data.
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
from shapely.geometry import LineString
from skimage import measure
import functools

from wavebreaking.utils.data_utils import get_dimension_attributes, check_argument_types
from wavebreaking.utils.index_utils import (
    iterate_time_dimension,
    iterate_contour_levels,
)


@check_argument_types(["data"], [xr.DataArray])
@get_dimension_attributes("data")
@iterate_time_dimension
@iterate_contour_levels
def calculate_contours(
    data, contour_levels, periodic_add=120, original_coordinates=True, *args, **kwargs
):
    """
    Calculate contour lines for a set of contour levels.
    The calculations are based on the measure.find_contours module.
    If periodic_add is provided, the data array is expanded in the longitudinal direction
    and undulations at the date border are captured correctly.
    Dimension names ("time_name", "lon_name", "lat_name"), size ("ntime", "nlon", "nlat")
    and resolution ("dlon", "dlat") can be passed as key=value argument.

    Parameters
    ----------
        data : xarray.DataArray
            data for the contour calculation
        contour_levels : array_like
            levels for contour calculation
        periodic_add: int or float, optional
            number of longitudes in degrees to expand the dataset to
            correctly capture undulations at the date border
            if the input field is not periodic, use periodic_add = 0
        original_coordinates: bool, optional
            if False, the array indices of the contour lines are returned
            instead of the original coordinates

    Returns
    -------
        contours: geopandas.GeoDataFrame
            GeoDataFrame containing different characteristics of the contours:
                * "date": date of the contour line
                * "level": level of the contour line
                * "closed": True if contour line is closed
                * "exp_lon": expansion in degrees of the contours in the longitudinal direction
                * "mean_lat": mean latitude of the contours
                * "geometry": LineString object with the contour coordinates in the format (x,y)
    """

    # select variable and time step for the contour calculation
    ds = data.sel({kwargs["time_name"]: kwargs["step"]})
    date = ds[kwargs["time_name"]].values

    # expand field for periodicity
    ds = xr.concat(
        [
            ds,
            ds.isel({kwargs["lon_name"]: slice(0, int(periodic_add / kwargs["dlon"]))}),
        ],
        dim=kwargs["lon_name"],
    )

    # get contours (indices in array coordinates)
    contours_from_measure = measure.find_contours(ds.values, kwargs["level"])

    # get indices and check closed and length
    contours_index_expanded = []
    closed = []
    for item in contours_from_measure:
        check_closed = all(item[0] == item[-1])
        indices = np.asarray(
            list(dict.fromkeys(map(tuple, np.round(item).astype("int"))))
        )[:, ::-1]
        check_len = len(indices) >= 4

        if check_len is True:
            contours_index_expanded.append(indices)
            closed.append(check_closed)

    def check_duplicates(list_of_arrays):
        """
        Check if there are cutoff duplicates due to the periodic expansion
        """
        temp = [
            np.c_[item[:, 0] % kwargs["nlon"], item[:, 1]] for item in list_of_arrays
        ]

        check = [
            (i1, i2)
            for (i1, e1), (i2, e2) in itertools.permutations(enumerate(temp), r=2)
            if set(map(tuple, e1)).issubset(set(map(tuple, e2)))
        ]
        drop = []
        lens = np.array([len(item) for item in temp])
        for indices in check:
            if lens[indices[0]] == lens[indices[1]]:
                drop.append(max(indices))
            else:
                drop.append(indices[np.argmin(lens[[indices[0], indices[1]]])])

        return list(set(drop))

    # check for duplicates
    drop = check_duplicates(contours_index_expanded)
    contours_index_expanded = [
        item for index, item in enumerate(contours_index_expanded) if index not in drop
    ]
    closed = [item for index, item in enumerate(closed) if index not in drop]

    def contour_to_dataframe(list_of_arrays):
        """
        Calculate different characteristics of the contour line.
        """

        exp_lon = [len(set(item[:, 0])) * kwargs["dlon"] for item in list_of_arrays]
        mean_lat = [np.round(item[:, 1].mean(), 2) for item in list_of_arrays]
        geo_mp = [LineString(coords) for coords in list_of_arrays]

        gdf = gpd.GeoDataFrame(
            pd.DataFrame(
                {
                    "date": date,
                    "level": kwargs["level"],
                    "closed": closed,
                    "exp_lon": exp_lon,
                    "mean_lat": mean_lat,
                }
            ),
            geometry=geo_mp,
        )

        return gdf

    if original_coordinates is False:
        # return contours in index coordinates as a geopandas.GeoDataFrame
        return contour_to_dataframe(contours_index_expanded)

    else:
        # select the original coordinates from the indices
        contours_coordinates_original = [
            np.c_[
                data[kwargs["lon_name"]].values[item[:, 0] % kwargs["nlon"]],
                data[kwargs["lat_name"]].values[item[:, 1]],
            ]
            for item in contours_index_expanded
        ]

        # drop duplicates
        contours_coordinates_original = [
            np.asarray(list(dict.fromkeys(map(tuple, item))))
            for item in contours_coordinates_original
        ]

        # return contours in original coordinates as a geopandas.GeoDataFrame
        return contour_to_dataframe(contours_coordinates_original)


def decorator_contour_calculation(func):
    """
    decorator to wrap the contour calculation around the index functions
    """

    @functools.wraps(func)
    def wrapper(data, contour_levels, periodic_add=120, *args, **kwargs):
        # pass contours to the decorated function as a key=value argument
        kwargs["contours"] = calculate_contours(
            data, contour_levels, periodic_add, original_coordinates=False
        )

        return func(data, contour_levels, periodic_add=periodic_add, *args, **kwargs)

    return wrapper
