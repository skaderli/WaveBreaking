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
import geopandas as gpd
from shapely.geometry import box
from tqdm import tqdm
import itertools as itertools
from sklearn.metrics import DistanceMetric

dist = DistanceMetric.get_metric("haversine")

from wavebreaking.utils.index_utils import properties_per_event
from wavebreaking.utils.data_utils import get_dimension_attributes, check_argument_types
from wavebreaking.indices.contour_index import decorator_contour_calculation


@check_argument_types(["data"], [xr.DataArray])
@get_dimension_attributes("data")
@decorator_contour_calculation
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
    Dimension names ("time_name", "lon_name", "lat_name"), size ("ntime", "nlon", "nlat")
    and resolution ("dlon", "dlat") can be passed as key=value argument.

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
                * "date": date of the overturning
                * "level": level of the contour line
                * "com": center of mass in the format (x,y)
                * "mean_var": mean of the variable used for the overturning calculations
                * "area": area of a overturning event;
                * "intensity": sum of the intensity (momentum flux)
                * "orientation": orientation of the most west- and eastward point
                * "geometry": (Multi)Polygon with coordinates in the format (x,y)
    """

    # filter contours from contour iteration
    contours = kwargs["contours"]
    contours = contours[contours.exp_lon == contours.exp_lon.max()].reset_index(
        drop=True
    )

    # loop over contours
    overturnings = []
    for index, series in tqdm(
        contours.iterrows(),
        total=contours.shape[0],
        desc="Calculating overturnings",
        leave=True,
        position=0,
    ):
        # calculate all overturning longitudes
        # (1) select the contour points and count the longitudes
        contour_index = pd.DataFrame(
            np.asarray(series.geometry.coords.xy).T, columns=["x", "y"]
        ).astype("int")
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

        # check if event is cyclonic by orientation
        def check_orientation(df):
            lat_west = data[kwargs["lat_name"]][
                contour_index[contour_index.x.eq(df.min_lon)].y.values[0]
            ].values
            lat_east = data[kwargs["lat_name"]][
                contour_index[contour_index.x.eq(df.max_lon)].y.values[-1]
            ].values

            if abs(lat_west) <= abs(lat_east):
                return "cyclonic"
            else:
                return "anticyclonic"

        orientation = pd.DataFrame(
            {"orientation": [check_orientation(row) for index, row in df_ot.iterrows()]}
        )
        # define Polygons
        dates = pd.DataFrame({"date": [series.date for i in range(0, len(df_ot))]})
        levels = pd.DataFrame({"level": [series.level for i in range(0, len(df_ot))]})
        polys = [
            box(row.min_lon, row.min_lat, row.max_lon, row.max_lat)
            for index, row in df_ot.iterrows()
        ]
        overturnings.append(
            gpd.GeoDataFrame(
                pd.concat([dates, levels, orientation], axis=1), geometry=polys
            )
        )

    # concat GeoDataFrames
    if len(overturnings) == 0:
        return gpd.GeoDataFrame()
    else:
        gdf = pd.concat(overturnings).reset_index(drop=True)

    # calculate properties and transform coordinates
    overturnings = [
        properties_per_event(data, row, intensity, periodic_add, *args, **kwargs)
        for (index, row) in tqdm(
            gdf.iterrows(),
            desc="Calculating properties  ",
            total=len(gdf),
            leave=True,
            position=0,
        )
    ]

    if len(overturnings) == 0:
        return gpd.GeoDataFrame()
    else:
        overturnings = pd.concat(overturnings).reset_index(drop=True)
        ot_gdf = gpd.GeoDataFrame(
            pd.concat([overturnings.iloc[:, :-1], gdf.orientation], axis=1),
            geometry=overturnings.geometry,
        )

        return ot_gdf
