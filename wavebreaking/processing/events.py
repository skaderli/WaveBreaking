""""""
"""
This file is part of WaveBreaking.

WaveBreaking provides indices to detect, classify
and track Rossby Wave Breaking (RWB) in climate and weather data.
The tool was developed during my master thesis at the University of Bern.
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---

Events post-processing functions
"""

__author__ = "Severin Kaderli"
__license__ = "MIT"
__email__ = "severin.kaderli@unibe.ch"

# import modules
import xarray as xr
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString, MultiPolygon
import shapely.vectorized
import numpy as np
from tqdm import tqdm
import itertools as itertools

from wavebreaking.utils.data_utils import (
    check_argument_types,
    check_empty_dataframes,
    get_dimension_attributes,
)
from wavebreaking.utils import index_utils


@check_argument_types(["data", "events"], [xr.DataArray, gpd.GeoDataFrame])
@check_empty_dataframes
@get_dimension_attributes("data")
def to_xarray(data, events, name="flag", *args, **kwargs):
    """
    Create xarray.DataArray from events stored in a geopandas.GeoDataFrame.
    Grid cells where an event is present are flagged with the value 1.
    Dimension names ("time_name", "lon_name", "lat_name"), size ("ntime", "nlon", "nlat")
    and resolution ("dlon", "dlat") can be passed as key=value argument.

    Parameters
    ----------
        data : xarray.DataArray
            data used for the index calculation
        events : geopandas.GeoDataFrame
            GeoDataFrame with the date and geometry for each event
        name : string, optional
            name of the xarray variable that is created

    Returns
    -------
        flag: xarray.DataArray
            Data with events flagged with the value 1
    """

    # get coordiantes of all events at the same time step
    def events_to_grid_points(geom, date):
        """
        Extract all grid cells that are enclosed by the path of a streamer
        """

        # get grid points
        x, y = np.meshgrid(data[kwargs["lon_name"]], data[kwargs["lat_name"]])
        x, y = x.flatten(), y.flatten()

        def polygon_to_grid(poly):
            # get interior points and border
            mask_contain = shapely.vectorized.contains(poly, x, y)
            mask_touch = shapely.vectorized.touches(poly, x, y)

            return np.c_[x, y][mask_contain | mask_touch]

        # get coords DataFrame
        if type(geom) == LineString:
            coords = np.asarray(geom.coords.xy).T
        elif type(geom) == MultiPolygon:
            coords = np.vstack([polygon_to_grid(poly) for poly in geom.geoms])
        else:
            coords = polygon_to_grid(geom)
        return pd.DataFrame(
            {"date": date, "x": coords[:, 0], "y": coords[:, 1]}
        ).drop_duplicates()

    # get all grid points of all events
    events_concat = pd.concat(
        [
            events_to_grid_points(row.geometry, row.date)
            for (index, row) in tqdm(
                events.iterrows(),
                desc="Converting to DataArray ",
                total=len(events),
                leave=True,
                position=0,
            )
        ]
    )

    # create empty xarray.Dataset with the same dimension as the original Dataset
    data_flagged = xr.zeros_like(data)

    # set coordinates of the events to the value 1
    data_flagged.loc[
        {
            kwargs["time_name"]: events_concat.date.to_xarray(),
            kwargs["lat_name"]: events_concat.y.to_xarray(),
            kwargs["lon_name"]: events_concat.x.to_xarray(),
        }
    ] = np.ones(len(events_concat))

    # change type and name
    data_flagged = data_flagged.astype("int8")
    data_flagged.name = name
    data_flagged.attrs["long_name"] = "flag wave breaking"

    return data_flagged


@check_argument_types(["events"], [pd.DataFrame])
@check_empty_dataframes
def track_events(events, time_range):
    """
    Temporal tracking of events.
    Events receive the same label if they spacially overlap at step t and at step t+1.

    Parameters
    ----------
        events : geopandas.GeoDataFrame
            GeoDataFrame with the date and coordinates of each identified event
        time_range: int or float
            Time range in hours for combining spatially overlapping events

    Returns
    -------
        events: geopandas.GeoDataFrame
            GeoDataFrame with label column showing the temporal connections
    """

    # resample event dates
    dates = [np.datetime64(date) for date in events.date]

    # reset index of events
    events = events.reset_index(drop=True)

    combine = []
    for date in dates:
        dates_dif = (dates - date).astype(int)
        check = (dates_dif >= 0) & (dates_dif <= time_range)
        index_combinations = itertools.combinations(events[check].index, r=2)
        combine.append(
            [
                combination
                for combination in index_combinations
                if events.iloc[combination[0]].geometry.intersects(
                    events.iloc[combination[1]].geometry
                )
            ]
        )

    # combine tracked indices to groups
    combine = index_utils.combine_shared(list(itertools.chain.from_iterable(combine)))

    # initiate label column
    events["label"] = events.index

    # assign lables to the events
    for item in combine:
        events.loc[item, "label"] = min(item)

    # select smallest possible index number for all events
    label = events.label.copy()
    for i in np.arange(len(set(events.label))):
        label[events.label == sorted(set(events.label))[i]] = i
    events.label = label

    # sort list by label and date and return geopandas.GeoDataFrame
    return events.sort_values(by=["label", "date"])
