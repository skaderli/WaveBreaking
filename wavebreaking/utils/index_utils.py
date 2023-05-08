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
from shapely.geometry import Polygon, MultiPolygon, GeometryCollection, box
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
            np.arange(0, kwargs["nlon"] + int(periodic_add / kwargs["dlon"])),
            np.arange(0, kwargs["nlat"]),
        )
        x, y = x.flatten(), y.flatten()

        # get interior points and border
        mask_contain = shapely.vectorized.contains(poly.buffer(0.5), x, y)

        return pd.DataFrame(np.c_[x, y][mask_contain], columns=["x", "y"]).astype("int")

    # get grid points
    points = polygon_to_grid_points(series.geometry)

    # select date and coordinates
    da = data.sel({kwargs["time_name"]: series.date}).values
    lons = data[kwargs["lon_name"]].values
    lats = data[kwargs["lat_name"]].values

    # Calculate cell weights following the definition by Daniel Steinfeld
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
    original_grid = box(0, 0, kwargs["nlon"] - 0.01, kwargs["nlat"])
    additional_grid = box(
        kwargs["nlon"],
        0,
        kwargs["nlon"] + int(periodic_add / kwargs["dlon"]),
        kwargs["nlat"],
    )

    def transform_coords(polygon):
        """
        Transform coordinates to original grid
        """
        coords = np.asarray(polygon.exterior.coords.xy).T.astype("int")
        coords = np.c_[coords[:, 0] % kwargs["nlon"], coords[:, 1]]
        return Polygon(np.c_[lons[coords[:, 0]], lats[coords[:, 1]]])

    if series.geometry.intersects(additional_grid):
        geoms_org = series.geometry.intersection(original_grid)
        geoms_add = series.geometry.intersection(additional_grid)

        def get_polygons(geoms):
            """
            Get polygons from intersections
            """
            if type(geoms) == Polygon:
                return [transform_coords(geoms)]
            elif type(geoms) == GeometryCollection:
                return [
                    transform_coords(geom)
                    for geom in geoms.geoms
                    if type(geom) == Polygon
                ]
            else:
                return []

        polys = get_polygons(geoms_org) + get_polygons(geoms_add)

        if len(polys) > 1:
            poly = MultiPolygon(polys)
        else:
            poly = polys[0]
    else:
        poly = transform_coords(series.geometry)

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
