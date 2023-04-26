"""
This file is part of WaveBreaking.

WaveBreaking provides indices to detect, classify
and track Rossby Wave Breaking (RWB) in climate and weather data.
The tool was developed during my master thesis at the University of Bern.
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---

Utility functions for the index calculation
"""

__author__ = "Severin Kaderli"
__license__ = "MIT"
__email__ = "severin.kaderli@unibe.ch"

# import modules
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
from shapely.geometry import MultiPoint
from tqdm import tqdm
import functools

# import logger
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)


def properties_per_timestep(date, level, lodf, data, intensity, periodic_add, *args, **kwargs):
    """
    Calculate several properties of the identified events.
    """

    # select time step
    ds = data.sel({kwargs["time_name"]: date})
    ds = ds.to_dataset(name="vars")

    # create variable "weight" representing the weighted area of each grid cell
    # this follows the weights definition used in
    # the ConTrack - Contour Tracking tool developed by Daniel Steinfeld
    weight_lat = np.cos(data[kwargs["lat_name"]].values * np.pi / 180)
    ds["weight"] = xr.DataArray(
        np.ones((kwargs["nlat"], kwargs["nlon"]))
        * np.array((111 * 1 * 111 * 1 * weight_lat)).astype(np.float32)[:, None],
        dims=[kwargs["lat_name"], kwargs["lon_name"]],
    )

    # add intensity if available
    if intensity is not None:
        ds["intensity"] = intensity.sel({kwargs["time_name"]: date})

    # expand for periodicity
    ds = xr.concat(
        [ds, ds.isel({kwargs["lon_name"]: slice(0, periodic_add)})],
        dim=kwargs["lon_name"],
    )

    def properties_per_event(df):
        """
        Calculate properties for each event.
        """

        # select grid cells of a specific event and calculate several properties
        temp = ds.isel(
            {kwargs["lat_name"]: df.y.to_xarray(), kwargs["lon_name"]: df.x.to_xarray()}
        )
        area = np.round(np.sum(temp.weight.values), 2)
        mean_var = np.round(np.average(temp.vars.values), 2)
        com = np.round(
            np.sum(df.multiply((temp.weight * temp.vars).values, axis=0), axis=0)
            / sum((temp.weight * temp.vars).values),
            2,
        )

        # transform coordinates to original coordinates
        x_original = data[kwargs["lon_name"]].values[
            df.x.values.astype("int") % kwargs["nlon"]
        ]
        y_original = data[kwargs["lat_name"]].values[
            df.y.values.astype("int") % kwargs["nlat"]
        ]
        geo_mp = MultiPoint(np.c_[x_original, y_original])

        com_x = data[kwargs["lon_name"]].values[
            (np.round(com.x, 2) % kwargs["nlon"]).astype("int")
        ]
        com_y = data[kwargs["lat_name"]].values[np.round(com.y, 2).astype("int")]

        prop_dict = {
            "date": date,
            "level": level,
            "com": (com_x, com_y),
            "mean_var": mean_var,
            "area": area,
        }

        if "intensity" in ds.variables:
            prop_dict["intensity"] = np.round(
                np.sum((temp.weight * temp.intensity).values) / area, 2
            )

        return [geo_mp, prop_dict]

    # store properties in a GeoDataFrame
    properties = [properties_per_event(item) for item in lodf]

    return gpd.GeoDataFrame(
        pd.DataFrame([item[1] for item in properties]),
        geometry=[item[0] for item in properties],
    )


def combine_shared(lst):
    """
    This is an internal function that combines all elements of a list
    that have at least one element in common.
    """

    elements = lst.copy()
    output = []
    while len(elements) > 0:
        first, *rest = elements
        first = set(first)

        lf = -1
        while len(first) > lf:
            lf = len(first)

            rest2 = []
            for r in rest:
                if len(first.intersection(set(r))) > 0:
                    first |= set(r)
                else:
                    rest2.append(r)
            rest = rest2

        output.append(list(first))
        elements = rest

    return output


def iterate_time_dimension(func):
    """
    decorator to iterate function over time dimension
    """

    @functools.wraps(func)
    def wrapper(data, contour_levels, *args, **kwargs):
        steps = data[kwargs["time_name"]]
        repeat_func = []

        for step in tqdm(steps, leave=True, position=0):
            kwargs["step"] = step
            repeat_func.append(func(data, contour_levels, *args, **kwargs))

        return pd.concat(repeat_func).reset_index(drop=True)

    return wrapper


def iterate_contour_levels(func):
    """
    decorator to iterate function over contour levels
    """

    @functools.wraps(func)
    def wrapper(data, contour_levels, *args, **kwargs):
        repeat_func = []

        try:
            iter(contour_levels)
        except Exception:
            contour_levels = [contour_levels]

        for level in contour_levels:
            kwargs["level"] = level
            repeat_func.append(func(data, contour_levels, *args, **kwargs))

        return pd.concat(repeat_func).reset_index(drop=True)

    return wrapper


def add_logger(message):
    """
    decorator to add a logger message
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            logger.info(message)

            return func(*args, **kwargs)

        return wrapper

    return decorator
