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
import numpy as np
import itertools as itertools
from sklearn.metrics import DistanceMetric

dist = DistanceMetric.get_metric("haversine")

from wavebreaking.utils.data_utils import (
    check_argument_types,
    check_empty_dataframes,
    get_dimension_attributes,
)
from wavebreaking.utils import index_utils


@check_argument_types(["data", "events"], [xr.DataArray, gpd.GeoDataFrame])
@check_empty_dataframes
@get_dimension_attributes("data")
def to_xarray(data, events, flag="ones", name="flag", *args, **kwargs):
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
        flag : string, optional
            column name of the events geopandas.GeoDataFrame
            flag is set where an event is present
            default value is "ones"
        name : string, optional
            name of the xarray variable that is created

    Returns
    -------
        flag: xarray.DataArray
            Data with events flagged with the value 1
    """

    # get grid points
    lon, lat = np.meshgrid(data[kwargs["lon_name"]], data[kwargs["lat_name"]])
    lonf, latf = lon.flatten(), lat.flatten()
    points = gpd.GeoDataFrame(
        pd.DataFrame({"lon": lonf, "lat": latf}),
        geometry=gpd.points_from_xy(lonf, latf),
    )

    # get coordinates of all events at the same time step
    buffer = events.copy()
    buffer.geometry = buffer.geometry.buffer(
        ((kwargs["dlon"] + kwargs["dlat"]) / 2) / 2
    )
    merged = gpd.sjoin(buffer, points, how="inner", predicate="contains").sort_index()

    # create empty xarray.Dataset with the same dimension as the original Dataset
    data_flagged = xr.zeros_like(data)

    # flag coordinates in DataSet
    if flag == "ones":
        set_val = np.ones(len(merged))
    else:
        try:
            set_val = merged[flag].values
        except KeyError:
            errmsg = "{} is not a column of the events geopandas.GeoDataFrame.".format(
                flag
            )
            raise KeyError(errmsg)

    data_flagged.loc[
        {
            kwargs["time_name"]: merged.date.to_xarray(),
            kwargs["lat_name"]: merged.lat.to_xarray(),
            kwargs["lon_name"]: merged.lon.to_xarray(),
        }
    ] = set_val

    # change type and name
    if flag == 'ones':
        data_flagged = data_flagged.astype("int8")
    data_flagged.name = name
    data_flagged.attrs["long_name"] = "flag wave breaking"

    return data_flagged


@check_argument_types(["events"], [pd.DataFrame])
@check_empty_dataframes
def track_events(
    events, time_range=None, method="by_overlap", buffer=0, overlap=0, distance=1000
):
    """
    Temporal tracking of events.
    Events receive the same label if they spatially overlap at step t
    and t + time_range.

    Parameters
    ----------
        events : geopandas.GeoDataFrame
            GeoDataFrame with the date and coordinates of each identified event
        time_range: int or float, optional
            Time range for temporally tracking the events. The units of
            time_range is hours if the type of the time dimension is np.datetime64.
            If not specified, the smallest time difference larger than zero is used.
        method : {"by_overlap", "by_distance"}, optional
            Method for temporally tracking the events:
                * "by_overlap": Events receive the same label if they spatially
                    overlap at step t and t + time_range.
                * "by_distance": Events receive the same label if their centre of mass
                    is closer than "distance"
        buffer : float, optional
            buffer around event polygon in degrees for the 'by_overlap' method
        overlap : float, optional
            minimum percentage of overlapping for the 'by_overlap' method
        distance : int or float, optional
            maximum distance in km between two events for the 'by_distance' method


    Returns
    -------
        events: geopandas.GeoDataFrame
            GeoDataFrame with label column showing the temporal coherence
    """

    # reset index of events
    events = events.reset_index(drop=True)

    # detect time range
    if time_range is None:
        date_dif = events.date.diff()
        time_range = date_dif[date_dif > pd.Timedelta(0)].min().total_seconds() / 3600

    # select events that are in range of time_range
    def get_range_combinations(events, index):
        """
        find events within the next steps that are in time range
        """
        if events.date.dtype == np.dtype("datetime64[ns]"):
            diffs = (events.date - events.date.iloc[index]).dt.total_seconds() / 3600
        else:
            diffs = abs(events.date - events.date.iloc[index])

        check = (diffs > 0) & (diffs <= time_range)

        return [(index, close) for close in events[check].index]

    range_comb = np.asarray(
        list(
            set(
                itertools.chain.from_iterable(
                    [get_range_combinations(events, index) for index in events.index]
                )
            )
        )
    )

    if len(range_comb) == 0:
        errmsg = "No events detected in the time range: {}".format(time_range)
        raise ValueError(errmsg)

    if method == "by_distance":
        # get centre of mass
        com1 = np.asarray(list(events.iloc[range_comb[:, 0]].com))
        com2 = np.asarray(list(events.iloc[range_comb[:, 1]].com))

        # calculate distance between coms
        dist_com = np.asarray(
            [dist.pairwise(np.radians([p1, p2]))[0, 1] for p1, p2 in zip(com1, com2)]
        )

        # check which coms are in range of 'distance'
        check_com = dist_com * 6371 < distance

        # select combinations
        combine = range_comb[check_com]

    elif method == "by_overlap":
        # select geometries that are in time range and add buffer
        geom1 = events.iloc[range_comb[:, 0]].geometry.buffer(buffer).make_valid()
        geom2 = events.iloc[range_comb[:, 1]].geometry.buffer(buffer).make_valid()

        # calculate and check the percentage of overlap
        inter = geom1.intersection(geom2, align=False)
        check_overlap = (
            inter.area.values
            / (geom2.area.values + geom1.area.values - inter.area.values)
            > overlap
        )

        # select combinations
        combine = range_comb[check_overlap]

    else:
        errmsg = "'{}' not supported as method!".format(method)
        hint = " Supported methods are 'by_overlap' and 'by_distance'"
        raise ValueError(errmsg + hint)

    # combine tracked indices to groups
    combine = index_utils.combine_shared(combine)

    # initiate label column
    events["label"] = events.index

    # assign labels to the events
    for item in combine:
        events.loc[item, "label"] = min(item)

    # select smallest possible index number for all events
    label = events.label.copy()
    for i in np.arange(len(set(events.label))):
        label[events.label == sorted(set(events.label))[i]] = i
    events.label = label

    # sort list by label and date and return geopandas.GeoDataFrame
    return events.sort_values(by=["label", "date"])
