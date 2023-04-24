"""
This file is part of WaveBreaking.

WaveBreaking provides indices to detect, classify
and track Rossby Wave Breaking (RWB) in climate and weather data.
The tool was developed during my master thesis at the University of Bern.
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---

Streamer calculation
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
from shapely.geometry import LineString, Polygon
import shapely.vectorized
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
@add_logger("Calculating streamers...")
def calculate_streamers(
    data,
    contour_level,
    geo_dis=800,
    cont_dis=1500,
    intensity=None,
    periodic_add=120,
    *args,
    **kwargs
):
    """
    Identify streamer structures based on the Streamer Index
    developed by Wernli and Sprenger (2007).
    The default parameters for the streamer identification
    are based on the study by Wernli and Sprenger (2007).

    Parameters
    ----------
        data : xarray.DataArray
            data for the contour and streamer calculation
        contour_level : int or float
            level of the contours
        geo_dis : int or float, optional
            Maximal geographic distance between two contour points that describe a streamer
        cont_dis : int or float, optional
            Minimal distance along the contour line between two points that describe a streamer
        intensity : xarray.DataArray, optional
            data for the intensity calculation (hint: use wb_spatial.calculate_momentum_flux)
        periodic_add: int or float, optional
            number of longitudes in degree to expand the dataset
            to correctly capture undulations at the date border
            if the input field is not periodic, use periodic_add = 0

    Returns
    -------
        streamers: geopandas.GeoDataFrame
            GeoDataFrame containing different characteristics of the streamers
                "date": date of the streamers
                "com": center of mass in the format (x,y)
                "mean_var": mean of the variable used for the streamer calculations
                "area": area of a streamer
                "intensity": sum of the intensity (momentum flux)
                over all streamer grid cells weighted with the corresponding area
                "geometry": MultiPoint object with the streamer coordinates in the format (x,y)
    """

    # filter contours
    contours = kwargs["contours"]
    contours = contours[contours.exp_lon == contours.exp_lon.max()].reset_index(
        drop=True
    )

    # loop over contours
    streamers = []
    for index, row in tqdm(
        contours.iterrows(), total=contours.shape[0], leave=True, position=0
    ):
        # calculate all possible basepoints combinations
        # (1) get coordinates of the contour points
        df_contour = pd.DataFrame([{"y": p.y, "x": p.x} for p in row.geometry.geoms])
        contour_line = LineString(row.geometry.geoms)

        # (2) calculate geographic distance between all contour coordinates
        geo_matrix = np.triu(
            (dist.pairwise((np.radians(df_contour.loc[:, ["y", "x"]]))) * 6371)
        )

        # (3) calculate the distance connecting all combinations of contour points
        on = np.insert(np.diagonal(geo_matrix, 1), 0, 0)
        on_mat = np.triu(np.tile(on, (len(on), 1)), k=1)
        cont_matrix = np.cumsum(on_mat, axis=1)

        # (4) get indices of all contour coordinates that fulfill both conditions
        check = (geo_matrix < geo_dis) * (cont_matrix > cont_dis)
        check = np.transpose(check.nonzero())

        # (5) store the coordinates of the point combinations in a DataFrame
        df1 = (
            df_contour[["y", "x"]]
            .iloc[check[:, 0]]
            .reset_index()
            .rename(columns={"y": "y1", "x": "x1", "index": "ind1"})
        )
        df2 = (
            df_contour[["y", "x"]]
            .iloc[check[:, 1]]
            .reset_index()
            .rename(columns={"y": "y2", "x": "x2", "index": "ind2"})
        )
        df_bp = pd.concat([df1, df2], axis=1)

        # (6) drop combinations that are invalid
        df_bp = df_bp.drop(df_bp[np.abs(df_bp.x1 - df_bp.x2) > 120].index)

        # apply several routines that neglect invalid basepoints
        def check_duplicates(df):
            """
            Check if there are basepoint duplicates
            due to the periodic expansion in the longitudinal direction
            """

            # drop duplicates after mapping the coordinates to the original grid
            temp = pd.concat(
                [df.y1, df.x1 % kwargs["nlat"], df.y2, df.x2 % kwargs["nlon"]], axis=1
            )
            temp = temp[temp.duplicated(keep=False)]

            if len(temp) == 0:
                check = []
            else:
                check = list(
                    np.asarray(
                        temp.groupby(list(temp))
                        .apply(lambda x: tuple(x.index))
                        .tolist()
                    )[:, 1]
                )

            return df.drop(check)

        def check_intersections(df):
            """
            Check for intersections of the basepoints with the contour line
            """

            # calculate line connecting each basepoint pair
            # and drop the ones that intersect with the contour line
            df["dline"] = [
                LineString([(row.x1, row.y1), (row.x2, row.y2)])
                for index, row in df.iterrows()
            ]
            check_intersections = [
                row.dline.touches(contour_line) for index, row in df.iterrows()
            ]

            return df[check_intersections]

        def check_overlapping(df):
            """
            Check for overlapping of the contour segments each described by basepoints
            """

            # calculate contour segment for each basepoint pair
            # and drop the pairs that are fully covered by another pair
            ranges = [range(int(r2.ind1), int(r2.ind2 + 1)) for i2, r2 in df.iterrows()]
            index_combinations = list(itertools.permutations(df.index, r=2))
            check_overlapping = set(
                [
                    item[0]
                    for item in index_combinations
                    if (
                        ranges[item[0]][0] in ranges[item[1]]
                        and ranges[item[0]][-1] in ranges[item[1]]
                    )
                ]
            )

            return df.drop(check_overlapping)

        def check_groups(df):
            """
            Check if there are still several basepoints describing the same streamer
            """

            # group basepoints that describe the same streamer
            index_combinations = np.asarray(
                list(itertools.combinations_with_replacement(df.index, r=2))
            )
            check_crossing = [
                df.iloc[item[0]].dline.intersects(df.iloc[item[1]].dline)
                for item in index_combinations
            ]

            groups = combine_shared(index_combinations[check_crossing])
            keep_index = [
                item[
                    np.argmax(
                        [
                            on[int(row.ind1): int(row.ind2) + 1].sum()
                            for index, row in df.iloc[item].iterrows()
                        ]
                    )
                ]
                for item in groups
            ]

            return df.iloc[keep_index]

        # apply all four routines and stop if no basepoints are left after a routine
        routines = [
            check_duplicates,
            check_intersections,
            check_overlapping,
            check_groups,
        ]

        i = 0
        while len(df_bp.index) > 1 and i <= 3:
            df_bp = routines[i](df_bp).reset_index(drop=True)
            i += 1

        def bp_to_grid(df):
            """
            Extract all grid cells that are enclosed by the path of a streamer
            """

            # map the streamer on the original grid
            x, y = np.meshgrid(
                np.arange(0, kwargs["nlon"] + periodic_add),
                np.arange(0, kwargs["nlat"]),
            )
            x, y = x.flatten(), y.flatten()
            mask = shapely.vectorized.contains(
                Polygon(np.array(df_contour[int(df.ind1): int(df.ind2) + 1])), y, x
            )

            return np.r_[
                df_contour[int(df.ind1): int(df.ind2) + 1].values, np.c_[y, x][mask]
            ]

        # return the result in a DataFrame
        if df_bp.empty:
            streamers.append(pd.DataFrame([]))
        else:
            streamer_grids = [
                pd.DataFrame(bp_to_grid(row), columns=["y", "x"]).astype("int")
                for index, row in df_bp.iterrows()
            ]
            streamers.append(
                properties_per_timestep(
                    row.date,
                    streamer_grids,
                    data,
                    intensity,
                    periodic_add,
                    *args,
                    **kwargs
                )
            )

    return pd.concat(streamers).reset_index(drop=True)
