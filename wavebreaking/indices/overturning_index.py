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
    contour_level,
    range_group=500,
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
        contour_level : int or float
            level of the contours
        range_group : int or float, optional
            Maximal distance in which two overturning events are grouped
        min_exp : int or float, optional
            Minimal longitudinal expansion of an overturning event
        intensity : xarray.DataArray, optional
            data for the intensity calculation (hint: use wb_spatial.calculate_momentum_flux)
        periodic_add: int or float, optional
            number of longitudes in degree to expand the dataset
            to correctly capture undulations at the date border
            if the input field is not periodic, use periodic_add = 0

    Returns
    -------
        overturnings: geopandas.GeoDataFrame
            GeoDataFrame containing different characteristics of the overturning events
                "date": date of the overturning
                "com": center of mass in the format (x,y)
                "mean_var": mean of the variable used for the overturning calculations
                "area": area of a overturning event
                "intensity": sum of the intensity (momentum flux)
                over all overturning grid cells weighted with the corresponding area
                "geometry": MultiPoint object with overturning coordinates in the format (x,y)
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
        df_contour = pd.DataFrame([{"y": p.y, "x": p.x} for p in row.geometry.geoms])
        lons, counts = np.unique(df_contour.x, return_counts=True)

        # (2) only keep longitudes that appear at least 3 times
        # and extract the first and last contour point at each of these longitudes
        ot_lons = lons[counts >= 3]
        ot_borders = [
            (
                df_contour[df_contour.x == lon].index[0],
                df_contour[df_contour.x == lon].index[-1],
            )
            for lon in ot_lons
        ]

        # (3) store date in a DataFrame
        df_ot = pd.DataFrame(
            [
                [
                    item[0],
                    df_contour.iloc[item[0]].y,
                    df_contour.iloc[item[0]].x,
                    item[1],
                    df_contour.iloc[item[1]].y,
                    df_contour.iloc[item[1]].x,
                ]
                for item in ot_borders
            ],
            columns=["ind1", "y1", "x1", "ind2", "y2", "x2"],
        )

        def check_duplicates(df):
            """
            Check if there are overturning longitude duplicates
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

        def check_joining(df):
            """
            Check for overlapping of the contour segment described by overturning borders
            """

            # combine overturning longitudes that share contour points in there segment
            ranges = [range(int(r2.ind1), int(r2.ind2) + 1) for i2, r2 in df.iterrows()]
            index_combinations = np.asarray(
                list(itertools.combinations_with_replacement(df.index, r=2))
            )
            check_join = [
                not elem
                for elem in [
                    set(ranges[item[0]]).isdisjoint(ranges[item[1]])
                    for item in index_combinations
                ]
            ]

            groups = combine_shared(index_combinations[check_join])
            idx = [
                (df.iloc[item].ind1.min(), df.iloc[item].ind2.max()) for item in groups
            ]

            return pd.DataFrame(
                [
                    [
                        item[0],
                        df_contour.iloc[item[0]].y,
                        df_contour.iloc[item[0]].x,
                        item[1],
                        df_contour.iloc[item[1]].y,
                        df_contour.iloc[item[1]].x,
                    ]
                    for item in idx
                ],
                columns=["ind1", "y1", "x1", "ind2", "y2", "x2"],
            )

        def check_distance(df):
            """
            Combine overturning events if they are closer than range_group
            """

            # calculate distance on the contour line between the overturning events
            geo_matrix = np.triu(
                (
                    dist.pairwise(
                        np.radians(df.loc[:, ["y2", "x2"]]),
                        np.radians(df.loc[:, ["y1", "x1"]]),
                    )
                    * 6371
                )
            )
            indices_check_dis = [
                (i, j)
                for i, j in itertools.combinations_with_replacement(df.index, r=2)
                if geo_matrix[i, j] < range_group
            ]

            groups = combine_shared(indices_check_dis + [(i, i) for i in df.index])
            idx = [
                (df.iloc[item].ind1.min(), df.iloc[item].ind2.max()) for item in groups
            ]

            return pd.DataFrame(
                [
                    [
                        item[0],
                        df_contour.iloc[item[0]].y,
                        df_contour.iloc[item[0]].x,
                        item[1],
                        df_contour.iloc[item[1]].y,
                        df_contour.iloc[item[1]].x,
                    ]
                    for item in idx
                ],
                columns=["ind1", "y1", "x1", "ind2", "y2", "x2"],
            )

        def check_expansion(df):
            """
            Drop overturning events that don't have a longitudinal expansion larger than min_exp
            """

            # calculate the longitudinal expansion
            check_exp = [
                (
                    max(df_contour[int(row.ind1): int(row.ind2) + 1].x)
                    - min(df_contour[int(row.ind1): int(row.ind2) + 1].x)
                )
                > min_exp
                for index, row in df.iterrows()
            ]
            return df[check_exp]

        # apply all four routines and stop if no overturning events are left after a routine
        routines = [check_duplicates, check_joining, check_distance, check_expansion]

        i = 0
        while len(df_ot.index) > 0 and i <= 3:
            df_ot = routines[i](df_ot).reset_index(drop=True)
            i += 1

        def ot_to_grid(df):
            """
            Get all grid cells that are enclosed by the rectangle describing an overturning event
            """

            # map the overturning events on the original grid,
            # represented by a rectangle around the event
            temp = df_contour[int(df.ind1): int(df.ind2) + 1]
            x, y = np.meshgrid(
                range(int(min(temp.x)), int(max(temp.x)) + 1),
                range(int(min(temp.y)), int(max(temp.y)) + 1),
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
                    overturning_grids,
                    data,
                    intensity,
                    periodic_add,
                    *args,
                    **kwargs
                )
            )

    return pd.concat(overturnings).reset_index(drop=True)
