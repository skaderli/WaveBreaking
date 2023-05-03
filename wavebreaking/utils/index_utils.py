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
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon
import shapely.vectorized
from tqdm import tqdm
import functools

# import logger
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)


def properties_per_event(data, series, intensity, periodic_add, *args, **kwargs):
    def polygon_to_grid_points(poly):
        """
        Extract all grid cells that are enclosed by the path of a polygon
        """
        # get grid points
        x, y = np.meshgrid(
            np.arange(0, kwargs["nlon"] + periodic_add),
            np.arange(0, kwargs["nlat"]),
        )
        x, y = x.flatten(), y.flatten()

        # get interior points and border
        mask_contain = shapely.vectorized.contains(poly.buffer(0.5), x, y)

        return pd.DataFrame(np.c_[x, y][mask_contain], columns=["x", "y"]).astype("int")

    # get grid points
    points = polygon_to_grid_points(series.geometry)

    # select date abd coordinates
    da = data.sel({kwargs["time_name"]: series.date}).values
    lons = data[kwargs["lon_name"]].values
    lats = data[kwargs["lat_name"]].values

    # Calculate cell weight following the definition by Daniel Steinfeld
    weight_lat = np.cos(data[kwargs["lat_name"]].values * np.pi / 180)
    weight = (
        np.ones((kwargs["nlat"], kwargs["nlon"]))
        * np.array((111 * 111 * weight_lat))[:, None]
    )

    # calculate properties
    areas = weight[points.y, points.x % kwargs["nlon"]]
    area = np.sum(areas)
    mean_var = np.dot(da[points.y, points.x % kwargs["nlon"]], areas) / area

    com = np.sum(points.multiply(areas, axis=0), axis=0) / area
    com = [(lons[int(com.x) % kwargs["nlon"]], lats[int(com.y)])]

    # transform coordinates and define polygons
    exterior = np.asarray(series.geometry.exterior.coords.xy).T.astype("int")

    if any(exterior[:, 0] >= kwargs["nlon"]):
        ext_east = exterior[exterior[:, 0] < kwargs["nlon"]]
        ext_west = exterior[exterior[:, 0] >= kwargs["nlon"]]

        if len(ext_east) < 4:
            poly = Polygon(
                np.c_[lons[ext_west[:, 0] % kwargs["nlon"]], lats[ext_west[:, 1]]]
            )
        elif len(ext_west) < 4:
            poly = Polygon(np.c_[lons[ext_east[:, 0]], lats[ext_east[:, 1]]])
        else:
            poly_east = Polygon(np.c_[lons[ext_east[:, 0]], lats[ext_east[:, 1]]])
            poly_west = Polygon(
                np.c_[lons[ext_west[:, 0] % kwargs["nlon"]], lats[ext_west[:, 1]]]
            )

            poly = MultiPolygon([poly_east, poly_west])

    else:
        poly = Polygon(np.c_[lons[exterior[:, 0]], lats[exterior[:, 1]]])

    # store properties in dict
    prop_dict = {
        "date": series.date,
        "level": series.level,
        "com": com,
        "mean_var": np.round(mean_var, 2),
        "area": np.round(area, 2),
    }

    # calculate intensity if provided
    if intensity is not None:
        intens = intensity.sel({kwargs["time_name"]: series.date}).values
        intens = np.dot(intens[points.y, points.x % kwargs["nlon"]], areas) / area
        prop_dict["intensity"] = np.round(intens, 2)

    # return GeoDataFrame
    return gpd.GeoDataFrame(
        pd.DataFrame(prop_dict),
        geometry=[poly],
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

        for step in tqdm(
            steps, desc="Calculating contours    ", leave=True, position=0
        ):
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
