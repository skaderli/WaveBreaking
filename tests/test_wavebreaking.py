"""
This file is part of WaveBreaking.

WaveBreaking provides indices to detect, classify
and track Rossby Wave Breaking (RWB) in climate and weather data.
The tool was developed during my master thesis at the University of Bern.
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---
Test classes for WaveBreaking
"""

__author__ = "Severin Kaderli"
__license__ = "MIT"
__email__ = "severin.kaderli@unibe.ch"

# import modules
import unittest
import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
import tqdm


# import and silent tqdm progress bar
def silent_tqdm(it, *args, **kwargs):
    return it


tqdm.tqdm = silent_tqdm


from wavebreaking.utils.data_utils import check_argument_types, get_dimension_attributes
from wavebreaking.utils.index_utils import (
    properties_per_event,
    iterate_time_dimension,
    iterate_contour_levels,
    combine_shared,
)
from wavebreaking.processing.spatial import (
    calculate_momentum_flux,
    calculate_smoothed_field,
)
from wavebreaking.processing.events import to_xarray, track_events
from wavebreaking.indices.contour_index import (
    calculate_contours,
    decorator_contour_calculation,
)
from wavebreaking.indices.streamer_index import calculate_streamers
from wavebreaking.indices.overturning_index import calculate_overturnings
from wavebreaking.indices.cutoff_index import calculate_cutoffs

data = xr.open_dataset("tests/data/demo_data.nc").isel(time=slice(0, 3))


class test_data_utils(unittest.TestCase):
    def test_check_argument_types(self):
        @check_argument_types(["data"], [xr.DataArray])
        def to_be_decorated(data, *args, **kwargs):
            return None

        self.assertRaises(TypeError, lambda: to_be_decorated(""))

    def test_get_dimension_attributes(self):
        @get_dimension_attributes("data")
        def to_be_decorated(data, *args, **kwargs):
            self.assertEqual(kwargs["time_name"], "time")
            self.assertEqual(kwargs["lon_name"], "lon")
            self.assertEqual(kwargs["lat_name"], "lat")

            self.assertEqual(kwargs["ntime"], 3)
            self.assertEqual(kwargs["nlon"], 360)
            self.assertEqual(kwargs["nlat"], 179)

            self.assertEqual(kwargs["dlon"], 1)
            self.assertEqual(kwargs["dlat"], 1)

        to_be_decorated(data.PV)


class test_index_utils(unittest.TestCase):
    def test_iterate_time_dimension(self):
        @get_dimension_attributes("data")
        @iterate_time_dimension
        def to_be_decorated(data, contour_levels, *args, **kwargs):
            return pd.DataFrame([kwargs["step"].values], columns=["time"])

        df_check = data.drop_vars("lev").time.to_dataframe().reset_index(drop=True)
        self.assertEqual(True, to_be_decorated(data.PV, 2).equals(df_check))

    def test_iterate_contour_levels(self):
        @get_dimension_attributes("data")
        @iterate_contour_levels
        def to_be_decorated(data, contour_levels, *args, **kwargs):
            return pd.DataFrame([kwargs["level"]], columns=["levels"])

        df_check = pd.DataFrame([-2, 2], columns=["levels"])
        self.assertEqual(True, to_be_decorated(data.PV, [-2, 2]).equals(df_check))

    def test_combine_shared(self):
        list_in = [[1, 2, 3], [2, 3, 4], [5, 6]]
        list_out = [[1, 2, 3, 4], [5, 6]]
        self.assertEqual(combine_shared(list_in), list_out)

    def test_properties_per_event(self):
        poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        date = "1959-06-03T12"
        level = 2
        gdf = gpd.GeoDataFrame(
            pd.DataFrame({"date": [date], "level": [level]}), geometry=[poly]
        )

        df_check = properties_per_event(
            data=data.PV,
            series=gdf.iloc[0],
            intensity=data.PV,
            periodic_add=120,
            time_name="time",
            lat_name="lat",
            nlat=179,
            nlon=360,
            lon_name="lon",
            dlon=1,
        )

        self.assertIs(type(df_check), gpd.GeoDataFrame)
        self.assertEqual(len(df_check.columns), 7)


class test_spatial(unittest.TestCase):
    def test_calculate_momentum_flux(self):
        self.assertIs(type(calculate_momentum_flux(data.U, data.V)), xr.DataArray)

    def test_calculate_smoothed_filed(self):
        self.assertIs(type(calculate_smoothed_field(data.PV, 10)), xr.DataArray)


class test_events(unittest.TestCase):
    def test_to_xarray(self):
        date = np.datetime64("1959-06-03T12")
        name = "test_flag"
        events = gpd.GeoDataFrame(
            pd.DataFrame([{"date": date}]),
            geometry=[Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])],
        )

        flag_data = to_xarray(data=data.PV, events=events, name=name)

        self.assertIs(type(flag_data), xr.DataArray)
        self.assertEqual(flag_data.sel(lon=5, lat=5, time=date).values, 1)
        self.assertEqual(flag_data.name, name)

    def test_track_events(self):
        date = np.datetime64("1959-06-03T12")
        events = gpd.GeoDataFrame(
            pd.DataFrame([{"date": date}]),
            geometry=[Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])],
        )

        tracked = track_events(events=events, time_range=24)

        self.assertIs(type(tracked), gpd.GeoDataFrame)
        self.assertEqual(tracked.iloc[0].label, 0)


class test_indices(unittest.TestCase):
    def test_contour_index(self):
        contours_coords = calculate_contours(
            data=data.PV,
            contour_levels=2,
            periodic_add=120,
            original_coordinates=True,
        )
        contours_coords = contours_coords[contours_coords.exp_lon == 360]

        contours_index = calculate_contours(
            data=data.PV,
            contour_levels=2,
            periodic_add=120,
            original_coordinates=False,
        )
        contours_index = contours_index[contours_index.exp_lon == 480]

        self.assertIs(type(contours_coords), gpd.GeoDataFrame)
        self.assertIs(type(contours_index), gpd.GeoDataFrame)

        self.assertEqual(
            min(np.asarray(contours_coords.iloc[0].geometry.coords.xy).T[:, 0]), -180
        )
        self.assertEqual(
            min(np.asarray(contours_index.iloc[0].geometry.coords.xy).T[:, 0]), 0
        )

        cols = ["date", "level", "closed", "exp_lon", "mean_lat", "geometry"]
        self.assertEqual(contours_coords.columns.to_list(), cols)
        self.assertEqual(contours_index.columns.to_list(), cols)

    def test_decorator_contour_calculation(self):
        @decorator_contour_calculation
        def to_be_decorated(*args, **kwargs):
            contours = kwargs["contours"]
            contours = contours[contours.exp_lon == 480]

            self.assertIs(type(contours), gpd.GeoDataFrame)
            self.assertEqual(
                min(np.asarray(contours.iloc[0].geometry.coords.xy).T[:, 0]), 0
            )
            self.assertEqual(
                contours.columns.to_list(),
                ["date", "level", "closed", "exp_lon", "mean_lat", "geometry"],
            )

        to_be_decorated(data.PV, contour_levels=2)

    def test_streamer_index(self):
        streamers = calculate_streamers(
            data=data.PV,
            contour_levels=2,
            geo_dis=800,
            cont_dis=1500,
            intensity=data.PV,
            periodic_add=120,
        )

        self.assertIs(type(streamers), gpd.GeoDataFrame)
        self.assertEqual(len(streamers), 39)
        self.assertEqual(
            streamers.columns.to_list(),
            ["date", "level", "com", "mean_var", "area", "intensity", "geometry"],
        )

    def test_overturning_index(self):
        overturnings = calculate_overturnings(
            data=data.PV,
            contour_levels=2,
            range_group=5,
            min_exp=5,
            intensity=data.PV,
            periodic_add=120,
        )

        self.assertIs(type(overturnings), gpd.GeoDataFrame)
        self.assertEqual(len(overturnings), 9)
        self.assertEqual(
            overturnings.columns.to_list(),
            [
                "date",
                "level",
                "com",
                "mean_var",
                "area",
                "intensity",
                "orientation",
                "geometry",
            ],
        )

    def test_cutoff_index(self):
        cutoffs = calculate_cutoffs(
            data=data.PV,
            contour_levels=2,
            min_exp=5,
            intensity=data.PV,
            periodic_add=120,
        )

        self.assertIs(type(cutoffs), gpd.GeoDataFrame)
        self.assertEqual(len(cutoffs), 20)
        self.assertEqual(
            cutoffs.columns.to_list(),
            ["date", "level", "com", "mean_var", "area", "intensity", "geometry"],
        )


# execute Test
if __name__ == "__main__":
    unittest.main()
