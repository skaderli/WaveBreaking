"""
This file is part of WaveBreaking. 

WaveBreaking provides indices to detect, classify and track Rossby Wave Breaking (RWB) in climate and weather data. 
The tool was developed during my master thesis at the University of Bern.
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---

Cutoff calculation
"""

__author__ = "Severin Kaderli"
__license__ = "MIT"
__email__ = "severin.kaderli@unibe.ch"

# import modules
import numpy as np
import pandas as pd
import itertools as itertools
from tqdm import tqdm
from shapely.geometry import Polygon
import shapely.vectorized
from sklearn.metrics import DistanceMetric
dist = DistanceMetric.get_metric('haversine')

from wavebreaking.utils.index_utils import *
from wavebreaking.utils.data_utils import *
from wavebreaking.indices.contour_index import decorator_contour_calculation

@get_dimension_attributes("data")
@decorator_contour_calculation
@add_logger("Calculating cutoffs...")
def calculate_cutoffs(data, contour_level, min_exp = 5, intensity = None, periodic_add = 120, *args, **kwargs):
    
    """
    Identify cutoff structures. 

    Parameters
    ----------    
        data : xarray.DataArray
            data for the contour and cutoff calculation
        contour_level : int or float
            level of the contours
        min_exp : int or float, optional
            Minimal longitudinal expansion of a cutoff event
        intensity : xarray.DataArray, optional
            data for the intensity calculation (hint: use wb_spatial.calculate_momentum_flux)
        periodic_add: int or float, optional
            number of longitudes in degree to expand the dataset to correctly capture undulations at the date border
            if the input field is not periodic, use periodic_add = 0

    Returns
    -------
        cutoffs: geopandas.GeoDataFrame
            GeoDataFrame containing different characteristics of the cutoffs events
                "date": date of the cutoffs
                "com": center of mass in the format (x,y)
                "mean_var": mean of the variable used for the contour calculation over all cutoff grid cells
                "area": area of a cutoff event
                "intensity": sum of the intensity (momentum flux) over all cutoff grid cells weighted with the corresponding area
            "geometry": MultiPoint object with the coordinate points of the cutoffs in the format (x,y)
    """
    
    #filter contours
    contours = kwargs["contours"]
    contours = contours[(contours.exp_lon < contours.exp_lon.max()) & (contours.exp_lon >= min_exp)].reset_index(drop=True)
    
    #loop over dates
    cutoffs = []
    for (name, group) in tqdm(contours.groupby("date"), leave = True, position = 0):

        def check_duplicates(df):

            """ 
            Check if there are cutoff duplicates due to the periodic expansion in the longitudinal direction 
            """

            temp = [[(p.y, p.x%kwargs["nlon"]) for p in row.geometry.geoms] for index, row in group.iterrows()]

            index_combinations = list(itertools.permutations(df.reset_index(drop=True).index, r=2))
            check = [item for item in index_combinations if set(temp[item[0]]).issubset(set(temp[item[1]]))]

            drop = []
            for item in check:
                lens = [len(temp[i]) for i in item]

                if lens[0] == lens[1]:
                    drop.append(max(item))
                else:
                    drop.append(item[np.argmin(lens)])

            return df[~df.reset_index(drop=True).index.isin(drop)].reset_index(drop=True)
        
        #check for duplicates
        df_co = check_duplicates(group)

        #drop cutoffs with less than 4 coordinates 
        df_co = df_co[[len(row.geometry.geoms)>=4 for index,row in df_co.iterrows()]]

        #drop cutoffs that are not closed
        df_co = df_co[[dist.pairwise((np.radians(np.asarray([(row.geometry.geoms[0].y, row.geometry.geoms[0].x), 
                                               (row.geometry.geoms[-1].y, row.geometry.geoms[-1].x)]))))[0,1]*6371 < 200 for index, row in df_co.iterrows()]]
        
        def cutoff_to_grid(df):

            """
            Extract all grid cells that are enclosed by the path of a cutoff
            """

            #map the cutoffs on the original grid
            x, y = np.meshgrid(np.arange(0,kwargs["nlon"]+periodic_add), np.arange(0,kwargs["nlat"]))
            x, y = x.flatten(), y.flatten()
            mask = shapely.vectorized.contains(Polygon(df.geometry.geoms),x,y)

            return np.r_[np.asarray([(p.x, p.y) for p in df.geometry.geoms]), np.c_[x,y][mask]]
        
        #return the result in a DataFrame  
        if df_co.empty:
            cutoffs.append(pd.DataFrame([]))
        else:
            cutoff_grids = [pd.DataFrame(cutoff_to_grid(row), columns = ["x", "y"]).astype("int") for index,row in df_co.iterrows()]
            cutoffs.append(properties_per_timestep(group.date.iloc[0],cutoff_grids, data, intensity, periodic_add, *args, **kwargs))

    return pd.concat(cutoffs).reset_index(drop=True)