""""""
"""
This file is part of WaveBreaking.

WaveBreaking provides indices to detect, classify
and track Rossby Wave Breaking (RWB) in climate and weather data.
The tool was developed during my master thesis at the University of Bern.
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---

Overturning calculation
"""

__author__ = "Severin Kaderli"
__license__ = "MIT"
__email__ = "severin.kaderli@unibe.ch"

# import modules
import numpy as np
import xarray as xr
import pandas as pd
from tqdm import tqdm
import itertools as itertools
from sklearn.metrics import DistanceMetric

dist = DistanceMetric.get_metric("haversine")

from wavebreaking.utils.index_utils import (
    properties_per_timestep,
    combine_shared,
    add_logger,
)
from wavebreaking.utils.data_utils import get_dimension_attributes, check_argument_types
from wavebreaking.indices.contour_index import decorator_contour_calculation


@check_argument_types(["data"], [xr.DataArray])
@get_dimension_attributes("data")
@decorator_contour_calculation
@add_logger("Calculating overturnings...")
def calculate_overturnings(
    data,
    contour_levels,
    range_group=5,
    min_exp=5,
    intensity=None,
    periodic_add=120,
    *args,
    **kwargs
):
    """
    Identify overturning structures based on the Overturning Index
    developed by Barnes and Hartmann (2012).
    The default parameters for the overturning identification
    are based on the study by Barnes and Hartmann (2012).
    The overturning region is represented by a rectangle
    enclosing all overturning contour points.

    Parameters
    ----------
        data : xarray.DataArray
            data for the contour and overturning calculation
        contour_levels : array_like
            levels for contour calculation
        range_group : int or float, optional
            Maximal degrees in the longitudinal direction in which two overturning are grouped
        min_exp : int or float, optional
            Minimal longitudinal expansion of an overturning event
        intensity : xarray.DataArray, optional
            data for the intensity calculation (hint: use wb_spatial.calculate_momentum_flux)
        periodic_add: int or float, optional
            number of longitudes in degrees to expand the dataset
            to correctly capture undulations at the date border
            if the input field is not periodic, use periodic_add = 0

    Returns
    -------
        overturnings: geopandas.GeoDataFrame
            GeoDataFrame containing different characteristics of the overturning events:
                "date": date of the overturning;
                "com": center of mass in the format (x,y);
                "mean_var": mean of the variable used for the overturning calculations;
                "area": area of a overturning event;
                "intensity": sum of the intensity (momentum flux);
                over all overturning grid cells weighted with the corresponding area;
                "geometry": MultiPoint object with overturning coordinates in the format (x,y).
    """

    # filter contours
    contours = kwargs["contours"]
    contours = contours[contours.exp_lon == contours.exp_lon.max()].reset_index(
        drop=True
    )

    # loop over contours
    overturnings = []
    for index, row in tqdm(
        contours.iterrows(), total=contours.shape[0], leave=True, position=0
    ):
        # calculate all overturning longitudes
        # (1) select the contour points and count the longitudes
        contour_index = pd.DataFrame([{"y": p.y, "x": p.x} for p in row.geometry.geoms])
        lons, counts = np.unique(contour_index.x, return_counts=True)

        # (2) only keep longitudes that appear at least 3 times
        # and create add same label if they are closer than
        # the parameter range_group
        ot_lons = pd.DataFrame({"lon": lons[counts >= 3]})
        ot_lons["label"] = (ot_lons.diff() > range_group / kwargs["dlon"]).cumsum()

        # (3) extract min and max longitude of each group
        groups = ot_lons.groupby("label")
        df_ot = groups.agg(["min", "max"]).astype("int").reset_index(drop=True)
        df_ot.columns = ["min_lon", "max_lon"]

        def check_duplicates(df):
            """
            Check if there are cutoff duplicates due to the periodic expansion
            """

            temp = [
                np.array(range(row.min_lon, row.max_lon + 1)) % kwargs["nlon"]
                for index, row in df.iterrows()
            ]

            index_combinations = list(itertools.permutations(df.index, r=2))
            check = [
                item
                for item in index_combinations
                if set(temp[item[0]]).issubset(set(temp[item[1]]))
            ]

            drop = []
            for item in check:
                lens = [len(temp[i]) for i in item]

                if lens[0] == lens[1]:
                    drop.append(max(item))
                else:
                    drop.append(item[np.argmin(lens)])

            return df[~df.reset_index(drop=True).index.isin(drop)]

        def check_expansion(df):
            exp_lon = df.max_lon - df.min_lon
            return df[exp_lon >= min_exp / kwargs["dlon"]]

        def find_lat_expansion(df):
            lats = [
                contour_index[
                    contour_index.x.isin(range(row.min_lon, row.max_lon + 1))
                ].y
                for index, row in df.iterrows()
            ]
            ot_lats = pd.DataFrame(
                [(item.min(), item.max()) for item in lats],
                columns=["min_lat", "max_lat"],
            ).astype("int")

            return pd.concat([df, ot_lats], axis=1)

        # apply all three routines and stop if no overturning events are left after a routine
        routines = [check_duplicates, check_expansion, find_lat_expansion]

        i = 0
        while len(df_ot.index) > 0 and i <= 2:
            df_ot = routines[i](df_ot).reset_index(drop=True)
            i += 1

        def ot_to_grid(df):
            """
            Get all grid cells that are enclosed by the rectangle describing an overturning event
            """

            # map the overturning events on the original grid,
            # represented by a rectangle around the event
            x, y = np.meshgrid(
                range(df.min_lon, df.max_lon + 1),
                range(df.min_lat, df.max_lat + 1),
            )
            x, y = x.flatten(), y.flatten()
            return np.vstack((y, x)).T

        # return the result in a DataFrame
        if df_ot.empty:
            overturnings.append(pd.DataFrame([]))
        else:
            overturning_grids = [
                pd.DataFrame(ot_to_grid(row), columns=["y", "x"])
                for index, row in df_ot.iterrows()
            ]
            overturnings.append(
                properties_per_timestep(
                    row.date,
                    row.level,
                    overturning_grids,
                    data,
                    intensity,
                    periodic_add,
                    *args,
                    **kwargs
                )
            )

    return pd.concat(overturnings).reset_index(drop=True)
