""""""
"""
This file is part of WaveBreaking.

WaveBreaking provides indices to detect, classify
and track Rossby Wave Breaking (RWB) in climate and weather data.
The tool was developed during my master thesis at the University of Bern.
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---

Cutoff calculation
"""

__author__ = "Severin Kaderli"
__license__ = "MIT"
__email__ = "severin.kaderli@unibe.ch"

# import modules
import xarray as xr
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
from shapely.geometry import Polygon

from wavebreaking.utils.index_utils import properties_per_event
from wavebreaking.utils.data_utils import get_dimension_attributes, check_argument_types
from wavebreaking.indices.contour_index import decorator_contour_calculation


@check_argument_types(["data"], [xr.DataArray])
@get_dimension_attributes("data")
@decorator_contour_calculation
def calculate_cutoffs(
    data, contour_level, min_exp=5, intensity=None, periodic_add=120, *args, **kwargs
):
    """
    Identify cutoff structures.
    Dimension names ("time_name", "lon_name", "lat_name"), size ("ntime", "nlon", "nlat")
    and resolution ("dlon", "dlat") can be passed as key=value argument.

    Parameters
    ----------
        data : xarray.DataArray
            data for the contour and cutoff calculation
        contour_levels : array_like
            levels for contour calculation
        min_exp : int or float, optional
            Minimal longitudinal expansion of a cutoff event
        intensity : xarray.DataArray, optional
            data for the intensity calculation (hint: use wb_spatial.calculate_momentum_flux)
        periodic_add: int or float, optional
            number of longitudes in degrees to expand the dataset
            to correctly capture undulations at the date border
            if the input field is not periodic, use periodic_add = 0

    Returns
    -------
        cutoffs: geopandas.GeoDataFrame
            GeoDataFrame containing different characteristics of the cutoffs events:
                * "date": date of the cutoffs
                * "level": level of the contour line
                * "com": center of mass in the format (x,y)
                * "mean_var": mean of the variable used for the contour calculation
                * "area": area of a cutoff event
                * "intensity": sum of the intensity (momentum flux)
                * "geometry": (Multi)Polygon with the coordinates in the format (x,y)
    """

    # filter contours from contour iteration
    contours = kwargs["contours"]
    contours = contours[
        (contours.exp_lon < contours.exp_lon.max())
        & (contours.exp_lon >= min_exp)
        & contours.closed
    ].reset_index(drop=True)

    # define Polygons
    polys = [Polygon(row.geometry) for index, row in contours.iterrows()]
    gdf = gpd.GeoDataFrame(contours[["date", "level"]], geometry=polys)

    # calculate properties and transform coordinates
    cutoffs = [
        properties_per_event(data, row, intensity, periodic_add, *args, **kwargs)
        for (index, row) in tqdm(
            gdf.iterrows(),
            desc="Calculating properties  ",
            total=len(gdf),
            leave=True,
            position=0,
        )
    ]

    if len(cutoffs) == 0:
        return gpd.GeoDataFrame()
    else:
        return pd.concat(cutoffs).reset_index(drop=True)
