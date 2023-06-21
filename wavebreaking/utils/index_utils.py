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
from shapely.ops import linemerge, unary_union, polygonize
from shapely.validation import make_valid
from tqdm import tqdm
import functools

# import logger
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)


def calculate_properties(events, data, intensity, periodic_add, **kwargs):
    # get all grid points
    x, y = np.meshgrid(
        np.arange(0, kwargs["nlon"] + int(periodic_add / kwargs["dlon"])),
        np.arange(0, kwargs["nlat"]),
    )
    x, y = x.flatten(), y.flatten()
    points = gpd.GeoDataFrame(
        pd.DataFrame({"x": x, "y": y}), geometry=gpd.points_from_xy(x, y)
    )

    # get coordinates of all events
    buffer = events.copy()
    buffer.geometry = buffer.geometry.buffer(
        ((kwargs["dlon"] + kwargs["dlat"]) / 2) / 2
    )
    merged = gpd.sjoin(buffer, points, how="inner", predicate="contains").sort_index()

    # get original coordinates
    merged["lon"] = data[kwargs["lon_name"]].values[merged.x % kwargs["nlon"]]
    merged["lat"] = data[kwargs["lat_name"]].values[merged.y]

    # calculate area equator in km^2
    area_cell = (
        np.round(6371 * 2 * np.pi / (360 / ((kwargs["dlon"] + kwargs["dlat"]) / 2)))
        ** 2
    )
    weight_lat = np.cos(np.radians(data[kwargs["lat_name"]].values)) * area_cell
    merged["areas"] = weight_lat[merged.y]

    # calculate mean_var, intensity and centre of mass
    merged["mean_var"] = (
        merged.areas
        * data.loc[
            {
                kwargs["time_name"]: merged.date.to_xarray(),
                kwargs["lat_name"]: merged.lat.to_xarray(),
                kwargs["lon_name"]: merged.lon.to_xarray(),
            }
        ]
    )

    if intensity is not None:
        merged["intensity"] = (
            merged.areas
            * intensity.loc[
                {
                    kwargs["time_name"]: merged.date.to_xarray(),
                    kwargs["lat_name"]: merged.lat.to_xarray(),
                    kwargs["lon_name"]: merged.lon.to_xarray(),
                }
            ]
        )
    else:
        merged["intensity"] = 0

    merged["x_com"] = merged.x * merged.areas
    merged["y_com"] = merged.y * merged.areas

    agg_merged = merged.groupby("id").agg(
        {
            "areas": "sum",
            "mean_var": "sum",
            "intensity": "sum",
            "x_com": "sum",
            "y_com": "sum",
        }
    )

    # calculate centre of mass
    com_x = data[kwargs["lon_name"]].values[
        (agg_merged.x_com / agg_merged.areas).astype("int") % kwargs["nlon"]
    ]
    com_y = data[kwargs["lat_name"]].values[
        (agg_merged.y_com / agg_merged.areas).astype("int")
    ]
    com = list(map(tuple, np.c_[com_x, com_y]))

    prop_dict = {
        "date": events.date,
        "level": events.level,
        "com": com,
        "mean_var": (agg_merged.mean_var / agg_merged.areas).round(2),
        "intensity": (agg_merged.intensity / agg_merged.areas).round(2),
        "event_area": agg_merged.areas.round(2),
    }

    # add orientation if available
    if "orientation" in events.columns:
        prop_dict["orientation"] = events.orientation

    return pd.DataFrame(prop_dict)


def transform_polygons(events, data, **kwargs):
    def transform_coords(polygon):
        """
        Transform coordinates to original grid
        """
        # get coordinates and check split and last meridian
        coords = np.asarray(polygon.exterior.coords.xy).T.astype("int")
        split = (coords[:, 0] >= kwargs["nlon"]).any()

        # transform coordinates
        coords = np.c_[coords[:, 0] % kwargs["nlon"], coords[:, 1]]

        return split, Polygon(
            np.c_[
                data[kwargs["lon_name"]][coords[:, 0]],
                data[kwargs["lat_name"]][coords[:, 1]],
            ]
        )

    def split_polys(polygon):
        """
        Split polygons at the last meridian
        """
        # define last meridian
        p00, p01 = [kwargs["nlon"] - 1, kwargs["nlat"]], [kwargs["nlon"] - 1, 0]
        p10, p11 = [kwargs["nlon"], 0], [kwargs["nlon"], kwargs["nlat"]]
        meridian = Polygon([p00, p01, p10, p11])

        # split polygons
        merged = linemerge([polygon.boundary, meridian.boundary])
        borders = unary_union(merged)
        polygons = [
            p
            for p in polygonize(borders)
            if not meridian.contains(p) and make_valid(polygon).contains(p)
        ]

        # transform if possible
        if len(polygons) == 0:
            return Polygon()
        elif len(polygons) == 1:
            return transform_coords(polygons[0])[1]
        else:
            polys = [transform_coords(p)[1] for p in polygons]
            return MultiPolygon(polys)

    # return GeoDataFrame
    gdf = gpd.GeoDataFrame(
        [transform_coords(row.geometry) for index, row in events.iterrows()],
        columns=["split", "geometry"],
    )
    gdf.loc[gdf.split, "geometry"] = [
        split_polys(row.geometry) for index, row in events[gdf.split].iterrows()
    ]

    return gdf


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
