#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Description
-----------
 
    This tool provides two indices to detect, classify, and track Rossby Wave Breaking (RWB) in weather and climate data. 
    The indices are based on the Streamer Index by Wernli and Sprenger (2007) and the Overturning Index by Barnes and Hartmann (2012).

Content
-------

    The following classes are available:
 
    wavebreaking: To create a wavebreaking object with functions to calculate wave breaking indices.
 
Examples
--------

    >>> filename = '/path/to/data.nc' 
    >>> wb = wavebreaking() 
    >>> wb.read_xarray(filename.nc)

Author
--------

    @author: severinkaderli (Severin Kaderli; severin.kaderli@unibe.ch)

"""

# =======
# import packages

# data
import numpy as np
import xarray as xr
import pandas as pd
from skimage import measure
import itertools as itertools
from tqdm import tqdm
import wrf as wrf
from shapely.geometry import Point, LineString, Polygon
import shapely.vectorized
from sklearn.metrics import DistanceMetric
dist = DistanceMetric.get_metric('haversine')

# plotting
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# logs
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# wavebreaking class
# =============================================================================

class wavebreaking(object):
    
    """
    wavebreaking class
    Author : Severin Kaderli, University of Bern , 2023
    Remark : The class structure is based on the ConTrack - Contour Tracking tool developed by Daniel Steinfeld
    """

    # number of instances initiated
    num_of_wavebreaking = 0

    def __init__(self, filename="", ds=None, **kwargs):
        
        """The constructor for wavebreaking class. Initialize a wavebreaking instance.
        
        If filename is given, try to load it directly.
        Arguments to the load function can be passed as key=value argument.
        
        Parameters
        ----------
            filename : string
                Datapath + filename
            ds : dataset
                xarray dataset
        """
        
        if not filename:
            if ds is None:
                self.ds = None
            else:
                self.ds = ds
            return
        
        try:
            self.ds = None
            self.read_xarray(filename, **kwargs)
        except Exception:
            try:
                self.ds = None
                self.read(filename, **kwargs)
            except:
                raise ValueError("Unkown fileformat. Known formats are netcdf.")
                
        wavebreaking.num_of_wavebreaking += 1
    
    def __repr__(self):
        try:
            string = "\
            Xarray dataset with {} time steps. \n\
            Available fields: {}".format(
                self.ntime, ", ".join(self.variables)
            )
        except AttributeError:
            # Assume it's an empty Blocking()
            string = "\
            Empty wavebreaking container.\n\
            Hint: use read() to load data."
        return string

    def __str__(self):
        return 'Class {}: \n{}'.format(self.__class__.__name__, self.ds)
  
    def __len__(self):
        return len(self.ds)
    
    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.ds, attr)

    
    def __getitem__(self, key):
        return self.ds[key]

    @property
    def ntime(self):
        """Return the number of time steps"""
        if len(self.ds.dims) != 3:
            logger.warning(
                "\nBe careful with the dimensions, "
                "you want dims = 3 and shape:\n"
                "(latitude, longitude, time)"
            )
            return self.ds.dims[self._get_name_time()]
        return self.ds.dims[self._get_name_time()]

    @property
    def variables(self):
        """Return the names of the variables"""
        return list(self.ds.data_vars)
    
    @property
    def dimensions(self):
        """Return the names of the dimensions"""
        return list(self.ds.dims)
    
    @property
    def grid(self):
        """Return the number of longitude and latitude grid"""
        if len(self.ds.dims) != 3:
            logger.warning(
                "\nBe careful with the dimensions, "
                "you want dims = 3 and shape:\n"
                "(latitude, longitude, time)"
            )
            return None
        string = "\
        latitude: {} \n\
        longitude: {}".format(
            self.ds.dims[self._get_name_latitude()], self.ds.dims[self._get_name_longitude()]
        ) 
        print(string)

    @property
    def dataset(self):
        """Return the dataset"""
        return self.ds

    
# ============================================================================
# Read / Import / Save data
# ============================================================================
    
    
    def read(self, filename, **kwargs):
        """
        Reads a file into a xarray dataset.
        
        Parameters
        ----------
            filename : string
                Valid path + filename
        """
        if self.ds is None:
            self.ds = xr.open_dataset(filename, **kwargs)
            logger.debug('read: {}'.format(self.__str__))
        else:
            errmsg = 'wavebreaking() is already set!'
            raise ValueError(errmsg)
            
    def read_xarray(self, ds):
        """
        Read an existing xarray data set.
        
        Parameter:
        ----------
            ds: data set
                Valid xarray data set.
        """
        if self.ds is None:
            if not isinstance(ds, xr.core.dataset.Dataset):
                errmsg = 'ds has to be a xarray data set!'
                raise ValueError(errmsg)
            self.ds = ds
            logger.debug('read_xarray: {}'.format(self.__str__))
        else:
            errmsg = 'wavebreaking() is already set!'
            raise ValueError(errmsg)

# ============================================================================
# Set up / Check dimensions
# ============================================================================

   
    def set_up(self,
               time_name=None,
               longitude_name=None,
               latitude_name=None,
               force=False,
    ):
        """
        Prepares the dataset for contour tracking. Does consistency checks
        and tests if all required information is available. Sets (automatically 
        or manually) internal variables and dimensions.
        Parameters
        ----------
            time_name : string, optional
                Name of time dimension. The default is None.
            longitude_name : string, optional
                Name of longitude dimension. The default is None.
            latitude_name : string, optional
                Name of latitude dimension. The default is None.
            force=False: bool, optional 
                Skip some consistency checks.
            write=True: bool, optional
                Print name of dimensions.
        Returns
        -------
            None.
        """
        
        # set dimensions
        if time_name is None:
            self._time_name = self._get_name_time()  
        else:
            self._time_name = time_name
        if longitude_name is None:
            self._longitude_name = self._get_name_longitude()
        else:
            self._longitude_name = longitude_name
        if latitude_name is None:
            self._latitude_name = self._get_name_latitude()
        else:
            self._latitude_name = latitude_name

        # set resolution
        if (self._longitude_name and self._latitude_name) is not None:
            self._dlon =  self._get_resolution(self._longitude_name, force=force)
            self._dlat =  self._get_resolution(self._latitude_name, force=force)

        if self._time_name is not None:
            self._dtime = self._get_resolution(self._time_name, force=force)

    
    def _get_name_time(self):
        """
        check for 'time' dimension and return name
        """
        # check unit
        for dim in self.ds.dims:
            if (('units' in self.ds[dim].attrs and
                'since' in self.ds[dim].attrs['units']) or 
                ('units' in self.ds[dim].encoding and
                 'since' in self.ds[dim].encoding['units']) or
                dim in ['time']):
                return dim
        # check dtype
        for dim in self.ds.variables:
            try:
                var = self.ds[dim].data[0]
            except IndexError:
                var = self.ds[dim].data
            if isinstance(var, np.datetime64):
                return dim   
        # no 'time' dimension found
        logger.warning(
            "\n 'time' dimension (dtype='datetime64[ns]') not found."
        )
        return None     


    def _get_name_longitude(self):
        """
        check for 'longitude' dimension and return name
        """
        for dim in self.ds.dims:
            if (('units' in self.ds[dim].attrs and
               self.ds[dim].attrs['units'] in ['degree_east', 'degrees_east']) or
               dim in ['lon', 'longitude', 'x']):
               return dim
        # no 'longitude' dimension found
        logger.warning(
            "\n 'longitude' dimension (unit='degrees_east') not found."
        )
        return None


    def _get_name_latitude(self):
        """
        check for 'latitude' dimension and return name
        """
        for dim in self.ds.dims:
            if (('units' in self.ds[dim].attrs  and
                self.ds[dim].attrs['units'] in ['degree_north', 'degrees_north']) or
                dim in ['lat', 'latitude', 'y']):
                return dim
        # no 'latitude' dimension found
        logger.warning(
            "\n 'latitude' dimension (unit='degrees_north') not found."
        )
        return None

    def _get_resolution(self, dim, force=False):
        """
        set spatial (lat/lon) and temporal (time) resolution
        """
        # time dimension in hours
        if dim == self._time_name:
            try:
                var = self.ds[dim].to_index()
                delta = np.unique((
                    self.ds[dim].to_index()[1:] - 
                    self.ds[dim].to_index()[:-1])
                    .astype('timedelta64[h]')
                )
            except AttributeError:  # dates outside of normal range
                # we can still move on if the unit is "days since ..."
                if ('units' in self.ds[dim].attrs and
                    'days' in self.ds[dim].attrs['units']):
                    var = self.ds[dim].data
                    delta = np.unique(var[1:] - var[:-1])
                else:
                    errmsg = 'Can not decode time with unit {}'.format(
                        self.ds[dim].attrs['units'])
                    raise ValueError(errmsg)
        # lat/lon dimension in Degree
        else:
            delta = abs(np.unique((
                self.ds[dim].data[1:] - 
                self.ds[dim].data[:-1])
            ))
        # check resolution
        if len(delta) > 1:
            errmsg = 'No regular grid found for dimension {}.\n\
            Hint: use set_up(force=True).'.format(dim)
            if force and dim != self._time_name:
                logging.warning(errmsg)
                logmsg = ' '.join(['force=True: using mean of non-equidistant',
                                   'grid {}'.format(delta)])
                logging.warning(logmsg)
                delta = round(delta.mean(), 2)
            else:
                if dim == self._time_name:
                    logging.warning(errmsg)
                else:
                    raise ValueError(errmsg)
        elif delta[0] == 0:
            errmsg = 'Two equivalent values found for dimension {}.'.format(
                dim)
            raise ValueError(errmsg)
        elif delta[0] < 0:
            errmsg = ' '.join(['{} not increasing. This should',
                                   'not happen?!']).format(dim)
            raise ValueError(errmsg)
            
        return delta

    
# ============================================================================
# Pre-processing functions
# ============================================================================


    def calculate_momentum_flux(self, variable_zonal, variable_meridional, dtime = "1D" ):
        
        """
        Calculate the momentum flux derived from the product of the deviations of both wind components from the zonal mean.
        
        Parameters
        ----------
            variable_zonal : string
                vairable name of the zonal wind variable.
            variable_meridional : string
                vairable name of the meridional wind variable.
            dtime : string, optional
                time period for the calculation of the zonal mean (e.g. "1D", "6H")
                
        Returns
        -------
            xarray.Dataset: float
                Dataset containing the variable "mflux"
        
        """

        if hasattr(self, '_time_name'):
            logger.info("Check dimensions... OK")
            pass
        else:
            logger.info("Set-up dimensions... DONE")
            self.set_up()

        logger.info('Calculating momentum flux...')

        #calculate zonal mean for each time step
        dzm_u = self.dataset[variable_zonal].resample({self._time_name:dtime}).mean(self._longitude_name)
        dzm_v = self.dataset[variable_meridional].resample({self._time_name:dtime}).mean(self._longitude_name)

        #calculate deviation from the zonal mean
        u_prime = self.dataset[variable_zonal]-dzm_u
        v_prime = self.dataset[variable_meridional]-dzm_v

        #mflux is given by the product of the deviation of both wind components
        self.dataset["mflux"] = u_prime*v_prime

        logger.info("Variable created: 'mflux'")
        
        
    def calculate_smoothed_field(self, variable, passes):
        
        """
        Calculate smoothed field based on a 5-point smoothing with double-weighted centre and multiple smoothing passes.
        The smoothing routine is based on the wrf.smooth2d function.
        
        Parameters
        ----------
            variable : string
                vairable name of the xarray.Dataset varialbe that should be smoothed
            passes : int or float
                number of smoothing passes of the 5-point smoothing
                
        Returns
        -------
            xarray.Dataset: float
                Dataset containing the the smoothed variable "smooth_variable"
        
        """
        
        if hasattr(self, '_time_name'):
            logger.info("Check dimensions... OK")
            pass
        else:
            logger.info("Set-up dimensions... DONE")
            self.set_up()

        logger.info('Calculating smoothed field...')

        #create wrap around the field the get a correct smoothing at the borders
        temp = self.dataset[variable].pad({self._longitude_name:5}, mode="wrap")
        #perform smoothing
        temp = wrf.smooth2d(temp,passes)
        #remove wrap from the smoothed field
        self.dataset[temp.name] = temp.isel({self._longitude_name:slice(5,-5)})
        
        #assign attributes
        self.dataset[temp.name] = self.dataset[temp.name].assign_attrs(self.dataset[variable].attrs)

        logger.info("Variable created: '{}'".format(temp.name))
    
    
# ============================================================================
# contour calculation
# ============================================================================

        
    def get_contours(self, variable, level, periodic_add = 120):
          
        """
        Calculate contour lines of a specific level at every time step. The calculations are based on the measure.find_contours functions.
        If periodic_add is provided, unclosed contours crossing the border in the longitudinal direction are closed. 
        
        Parameters
        ----------
            variable : string
                variable name of the xarray.Dataset whereof the contour lines should be calculated
            level : int or float
                level of the contour lines
            periodic_add: int or float, optional
                number of longitudes that are added to the dataset to correctly identify contours crossing the longitudinal border
                if the input field is not periodic, use periodic_add = 0
                
        Returns
        -------
            pd.DataFrame: DataFrame
                DataFrame providing different characteristics of the contours accessible at wavebreaking.contours
                    "date": date of the contour line in time steps of the original variable
                    "coordinates": coordinate points of the contour line (segement)
                    "level": level of the contour lines
                    "exp_lon": expansion of the contour line (segment) in the longitudinal direction
                    "mean_lat": mean latitude of the coordinates of a contour line (segment)
                    
                Information:
                The contours are mapped on the original grid of the variable. 
                The extended contours (if periodic_add > 0) used for the index calculation are accessible at wavebreaking._contours_wb_calc
        
        """
  
        if hasattr(self, '_time_name'):
            logger.info("Check dimensions... OK")
            pass
        else:
            logger.info("Set-up dimensions... DONE")
            self.set_up()
        
        logger.info('Calculating contour lines...')
        
        #set up progress bar
        times = self.dataset[self._time_name]
        progress = tqdm(range(0,len(times)), leave = True, position = 0)
        
        #store contour parameters
        self._periodic_add = periodic_add
        self._contour_level = level
        self._contour_variable = variable
        
        #calculate contours for each time step
        contours = [self.get_contour_timestep(step) for step,p in zip(times,progress)]
        
        #store contours that are mapped on the original grid
        self.contours = pd.concat([item[0] for item in contours]).reset_index(drop=True)
        
        #store contours that are used for the index calculation
        contours_wb_calc = pd.concat([item[1] for item in contours])
        self._contours_wb_calc = contours_wb_calc[contours_wb_calc.exp_lon == contours_wb_calc.exp_lon.max()].reset_index(drop=True)
        
        logger.info('{} contour(s) identified!'.format(len(self.contours)))
        
    def get_contour_timestep(self, step):
        """
            Calculate contour lines for one time step.
        """

        #select variable and time step for the contour calculation
        ds = self.dataset.sel({self._time_name:step})[self._contour_variable]
        #expand field for periodicity
        ds = xr.concat([ds, ds.isel({self._longitude_name:slice(0,self._periodic_add)})], dim=self._longitude_name)

        #get contours (indices in array coordinates)
        contours_from_measure = measure.find_contours(ds.values, self._contour_level)

        #contours for wave breaking calculation
        contours_index_expanded = [np.asarray(list(dict.fromkeys(map(tuple,np.round(item).astype("int"))))) for item in contours_from_measure]

        #map in added indices to the original indices of the grid
        contours_index_original = [np.c_[item[:,0]%self.dims[self._latitude_name], item[:,1]%self.dims[self._longitude_name]] for item in contours_index_expanded]
        contours_index_original =  [list(dict.fromkeys(map(tuple,item))) for item in contours_index_original]
        
        #drop duplicates that have been identified due to the periodicity
        index = np.arange(0,len(contours_index_original))
        index_combination = list(itertools.combinations(index, r=2))

        drop = []
        for combination in index_combination:
            if set(contours_index_original[combination[0]]) < set(contours_index_original[combination[1]]):
                drop.append(combination[0]) 
            if set(contours_index_original[combination[1]]) < set(contours_index_original[combination[0]]):
                drop.append(combination[1])
            if set(contours_index_original[combination[1]]) == set(contours_index_original[combination[0]]):
                drop.append(combination[1])

        contours_index_original = [np.asarray(contours_index_original[i]) for i in index if i not in drop]
        
        #select the original coordinates from the indices
        contour_coords_original = [np.c_[self.dataset[self._latitude_name].values[item[:,0].astype("int")], self.dataset[self._longitude_name].values[item[:,1].astype("int")]] for item in contours_index_original]
        
        #return contour in a DataFrame and calculate some characteristics
        def contour_to_df(list_of_arrays):
            """
                Calculate the different characteristics of the contour line.
            """
            date = pd.to_datetime(str(step.values)).strftime('%Y-%m-%dT%H')
            exp_lon = [len(set(item[:,1])) for item in list_of_arrays]
            mean_lat = [np.round(item[:,0].mean(),2) for item in list_of_arrays]
            coords = [list(map(tuple,item)) for item in list_of_arrays]

            return pd.DataFrame([{"date":date, "coordinates":c, "level":self._contour_level, "exp_lon":x, "mean_lat":m} for c,x,m in zip(coords, exp_lon, mean_lat) if len(c)>1])

        return contour_to_df(contour_coords_original), contour_to_df(contours_index_expanded)

        
# ============================================================================
# index calculations    
# ============================================================================
        
    
    def get_streamers(self, geo_dis = 800, cont_dis = 1500):
        
        """
        Identify streamer structures based on the Streamer Index developed by Wernli and Sprenger (2007). 
        The streamer calculation is based on the contour lines from wb.get_contours accessible at wavebreaking._contours_wb_calc.
        The default parameters for the streamer identification are based on the study by Wernli and Sprenger (2007).
        
        Parameters
        ----------
            geo_dis : int or float, optional
                Maximal geographic distance between two contour points that describe a streamer
            cont_dis : int or float, optional
                Minimal distance along the contour line between two contour points that describe a streamer
                
        Returns
        -------
            pd.DataFrame: DataFrame
                DataFrame providing different characteristics of the streamers accessible at wavebreaking.streamers
                    "date": date of the streamer in time steps of the original variable
                    "coordinates": coordinate points of the grid cells representing the streamer in the format (y,x)
                    "mean_var": mean of the variable used for the contour calculation over all streamer grid cell
                    "area": area of a streamer
                    "intensity": sum of the intensity (momentum flux) over all streamer grid cells weighted with the corresponding area
                    "com": center of mass of a streamer in the format (y,x)
                    
                Remark:
                The contours are mapped on the original grid of the variable. 
                The extended contours (if periodic_add > 0) used for the index calculation are accessible at wavebreaking._contours_wb_calc
        
        """
        
        if hasattr(self, '_contours_wb_calc'):
            logger.info("Check contours... OK")
            pass
        else:
            errmsg = "Check contours... No contours found! Calculate contours first using wavebreaking.get_contours()"
            raise ValueError(errmsg)
            
        if "mflux" in self.variables:
            logger.info("Momentum flux variable detected. Intensity will be calculated.")
        else:
            logger.info("No momentum flux detected. Intensity will not be calculated.")
       
        #create variable "weight" representing the weighted area of each grid cell
        #this follows the weights definition used in the ConTrack - Contour Tracking tool developed by Daniel Steinfeld
        weight_lat = np.cos(self.dataset[self._latitude_name].values*np.pi/180)
        self.dataset["weight"] = xr.DataArray(np.ones((self.dims[self._latitude_name], self.dims[self._longitude_name])) * np.array((111 * 1 * 111 * 1 * weight_lat)).astype(np.float32)[:, None], dims = [self._latitude_name, self._longitude_name])
        
        logger.info('Calculating streamers...')
        
        #set up progress bar
        times = self.dataset[self._time_name]
        progress = tqdm(range(0,len(times)), leave = True, position = 0)
        
        #calculate streamers for each time step
        self.streamers = pd.concat([self.get_streamers_timestep(row, geo_dis, cont_dis) for (index,row),p in zip(self._contours_wb_calc.iterrows(),progress)]).reset_index(drop=True)
        
        logger.info('{} streamer(s) identified!'.format(len(self.streamers)))
        
    def get_streamers_timestep(self, contour, geo_dis, cont_dis):
        """
            Calculate streamers for one time step. 
        """
        #calculate all possible basepoints combinations
        #(1) get coordinates of the contour points
        df_contour = pd.DataFrame(contour.coordinates, columns = ["y", "x"])
        contour_line = LineString([Point(row.x, row.y) for index, row in df_contour.iterrows()]) 

        #(2) calculate geographic distance between all contour coordinates
        geo_matrix = np.triu((dist.pairwise((np.radians(df_contour.loc[:,["y", "x"]])))*6371))

        #(3) calculate the distance connecting all combinations of contour points
        on = np.insert(np.diagonal(geo_matrix,1),0,0)
        on_mat = np.triu(np.tile(on,(len(on),1)), k=1)
        cont_matrix = np.cumsum(on_mat, axis = 1)

        #(4) get indices of all contour coordinates that fulfill both conditions
        check = (geo_matrix<geo_dis) * (cont_matrix>cont_dis)
        check = np.transpose(check.nonzero())

        #(5) store the coordinates of the point combinations in a DataFrame
        df1 = df_contour[["y","x"]].iloc[check[:,0]].reset_index().rename(columns = {"y":"y1", "x":"x1", "index":"ind1"})
        df2 = df_contour[["y","x"]].iloc[check[:,1]].reset_index().rename(columns = {"y":"y2", "x":"x2", "index":"ind2"})
        df_bp = pd.concat([df1, df2], axis = 1)
        
        #(6) drop combinations that are invalid
        df_bp = df_bp.drop(df_bp[np.abs(df_bp.x1-df_bp.x2)>120].index)

        #apply several routines that neglect invalid basepoints
        def check_duplicates(df):
            """
                Check if there are basepoint duplicates due to the periodic expansion in the longitudinal direction
            """
            #drop duplicates after mapping the coordinates to the original grid
            temp = pd.concat([df.y1, df.x1%self.dims[self._longitude_name], df.y2, df.x2%self.dims[self._longitude_name]], axis = 1)
            temp = temp[temp.duplicated(keep=False)]
            if len(temp) == 0:
                check = []
            else:
                check = list(np.asarray(temp.groupby(list(temp)).apply(lambda x: tuple(x.index)).tolist())[:,1])

            return df.drop(check)

        def check_intersections(df):
            """
                Check for intersections of the basepoints with the contour line 
            """
            #calculate line connecting each basepoint pair and drop the ones that intersect with the contour line
            df["dline"] = [LineString([(row.x1, row.y1),(row.x2, row.y2)]) for index, row in df.iterrows()]
            check_intersections = [row.dline.touches(contour_line) for index, row in df.iterrows()]

            return df[check_intersections]

        def check_overlapping(df):
            """
                Check for overlapping of the contour segments each described by basepoints
            """
            #calculate contour segment for each basepoint pair and drop the pairs that are fully covered by another pair
            ranges = [range(int(r2.ind1),int(r2.ind2+1)) for i2,r2 in df.iterrows()]
            index_combinations = list(itertools.permutations(df.index, r=2))
            check_overlapping = set([item[0] for item in index_combinations if (ranges[item[0]][0] in ranges[item[1]] and ranges[item[0]][-1] in ranges[item[1]])])

            return df.drop(check_overlapping)

        def check_groups(df):
            """
                Check if there are still several basepoints describing the same streamer
            """
            #group basepoints that describe the same streamer
            index_combinations = np.asarray(list(itertools.combinations_with_replacement(df.index, r=2)))
            check_crossing = [df.iloc[item[0]].dline.intersects(df.iloc[item[1]].dline) for item in index_combinations]

            groups = self.combine_shared(index_combinations[check_crossing])
            keep_index = [item[np.argmax([on[int(row.ind1):int(row.ind2)+1].sum() for index, row in df.iloc[item].iterrows() ])] for item in groups]

            return df.iloc[keep_index]

        #apply all four routines and stop if no basepoints are left after a routine
        routines = [check_duplicates, check_intersections, check_overlapping, check_groups]
        i = -1
        while len(df_bp.index) > 1:
            i+=1
            df_bp = routines[i](df_bp).reset_index(drop=True)
            if i == 3:
                break

        #return the result in a DataFrame  
        if df_bp.empty:
            return pd.DataFrame([])
        else:
            def bp_to_grid(df):
                """
                    Extract all grid cells that are enclosed by the path of a streamer
                """
                #map the streamer on the original grid
                x, y = np.meshgrid(np.arange(0,self.dims[self._longitude_name]+self._periodic_add), np.arange(0,self.dims[self._latitude_name]))
                x, y = x.flatten(), y.flatten()
                mask = shapely.vectorized.contains(Polygon(np.array(df_contour[int(df.ind1):int(df.ind2)+1])),y,x)
                
                return np.r_[df_contour[int(df.ind1):int(df.ind2)+1].values, np.c_[y,x][mask]]

            streamer_grids = [pd.DataFrame(bp_to_grid(row), columns = ["y", "x"])  for index,row in df_bp.iterrows()]
            return self.get_events(contour.date, streamer_grids)
    
    def get_overturnings(self, range_group = 500, min_exp = 5):
        
        """
        Identify overturning structures based on the Overturning Index developed by Barnes and Hartmann (2012). 
        The overturning calculation is based on the contour lines from wb.get_contours() accessible at wavebreaking._contours_wb_calc
        The default parameters for the overturning identification are based on the study by Barnes and Hartmann (2012).
        The overturning region is represented by a rectangle enclosing all overturing contour points
        
        Parameters
        ----------
            range_group : int or float, optional
                Maximal distance in which two overturning events are grouped
            min_exp : int or float, optional
                Minimal longitudinal expansion of an overturning event
                
        Returns
        -------
            pd.DataFrame: DataFrame
                DataFrame providing different characteristics of the overturning event accessible at wavebreaking.overturnings
                    "date": date of the overturning in time steps of the original variable
                    "coordinates": coordinate points of the grid cells representing the overturning event in the format (y,x)
                    "mean_var": mean of the variable used for the contour calculation over all streamer grid cell
                    "area": area of an overturning event
                    "intensity": sum of the intensity (momentum flux) over all overturning grid cells weighted with the corresponding area
                    "com": center of mass of an overturning event in the format (y,x)
                    
                Remark:
                The contours are mapped on the original grid of the variable. 
                The extended contours (if periodic_add > 0) used for the index calculation are accessible at wavebreaking._contours_wb_calc
        
        """
        
        if hasattr(self, '_contours_wb_calc'):
            logger.info("Check contours... OK")
            pass
        else:
            errmsg = "Check contours... No contours found! Calculate contours first using wavebreaking.get_contours()"
            raise ValueError(errmsg)
            
        if "mflux" in self.variables:
            logger.info("Momentum flux variable detected. Intensity will be calculated.")
        else:
            logger.info("No momentum flux detected. Intensity will not be calculated.")
       
        #create variable "weight" representing the weighted area of each grid cell
        #this follows the weights definition used in the ConTrack - Contour Tracking tool developed by Daniel Steinfeld
        weight_lat = np.cos(self.dataset[self._latitude_name].values*np.pi/180)
        self.dataset["weight"] = xr.DataArray(np.ones((self.dims[self._latitude_name], self.dims[self._longitude_name])) * np.array((111 * 1 * 111 * 1 * weight_lat)).astype(np.float32)[:, None], dims = [self._latitude_name, self._longitude_name])
        
        logger.info('Calculating overturning events...')
        
        #set up progress bar
        times = self.dataset[self._time_name]
        progress = tqdm(range(0,len(times)), leave = True, position = 0)
        
        #calculate overturning events for each time step
        self.overturnings = pd.concat([self.get_overturnings_timestep(row, range_group, min_exp) for (index,row),p in zip(self._contours_wb_calc.iterrows(),progress)]).reset_index(drop=True)
        
        logger.info('{} overturning event(s) identified!'.format(len(self.overturnings)))
    
    def get_overturnings_timestep(self, contour, range_group, min_exp):
        """
            Calculate overturning events for one time step. 
        """
        #calculate all overturning longitudes 
        #(1) select the contour points and count the longitudes
        df_contour = pd.DataFrame(contour.coordinates, columns = ["y", "x"])  
        lons, counts = np.unique(df_contour.x, return_counts = True)
        
        #(2) only keep longitudes that appear at least 3 times and extract the first and last contour point at each of these longitudes
        ot_lons = lons[counts >= 3]
        ot_borders = [(df_contour[df_contour.x == lon].index[0], df_contour[df_contour.x == lon].index[-1]) for lon in ot_lons]

        #(3) store date in a DataFrame
        df_ot =  pd.DataFrame([[item[0], df_contour.iloc[item[0]].y, df_contour.iloc[item[0]].x, 
                            item[1],  df_contour.iloc[item[1]].y, df_contour.iloc[item[1]].x] for item in ot_borders], columns = ["ind1","y1","x1","ind2","y2","x2"])

        #apply several routines to process the overturning events
        def check_duplicates(df):
            """
                Check if there are overturning longitude duplicates due to the periodic expansion in the longitudinal direction
            """
            #drop duplicates after mapping the coordinates to the original grid
            temp = pd.concat([df.y1, df.x1%self.dims[self._longitude_name], df.y2, df.x2%self.dims[self._longitude_name]], axis = 1)
            temp = temp[temp.duplicated(keep=False)]
            if len(temp) == 0:
                check = []
            else:
                check = list(np.asarray(temp.groupby(list(temp)).apply(lambda x: tuple(x.index)).tolist())[:,1])

            return df.drop(check)

        def check_joining(df):
            """
                Check for overlapping of the contour segment described by overturning borders
            """
            #combine overturning longitudes that share contour points in there segment
            ranges = [range(int(r2.ind1),int(r2.ind2)+1) for i2,r2 in df.iterrows()]
            index_combinations = np.asarray(list(itertools.combinations_with_replacement(df.index, r=2)))
            check_join = [not elem for elem in [set(ranges[item[0]]).isdisjoint(ranges[item[1]]) for item in index_combinations]]

            groups = self.combine_shared(index_combinations[check_join])
            idx = [(df_ot.iloc[item].ind1.min(), df_ot.iloc[item].ind2.max()) for item in groups]

            return pd.DataFrame([[item[0], df_contour.iloc[item[0]].y, df_contour.iloc[item[0]].x, 
                                  item[1], df_contour.iloc[item[1]].y, df_contour.iloc[item[1]].x] for item in idx], columns = ["ind1","y1","x1","ind2","y2","x2"])

        def check_distance(df):
            """
                Combine overturning events if they are closer than range_group
            """
            #calculate distance on the contour line between the overturning events
            geo_matrix = np.triu((dist.pairwise(np.radians(df.loc[:,['y2','x2']]), np.radians(df.loc[:,['y1','x1']]))*6371))
            indices_check_dis = [(i,j) for i,j in itertools.combinations_with_replacement(df.index, r=2) if geo_matrix[i,j] < range_group]

            groups = self.combine_shared(indices_check_dis+[(i,i) for i in df.index])
            idx = [(df.iloc[item].ind1.min(), df.iloc[item].ind2.max()) for item in groups]

            return pd.DataFrame([[item[0], df_contour.iloc[item[0]].y, df_contour.iloc[item[0]].x, 
                                  item[1], df_contour.iloc[item[1]].y, df_contour.iloc[item[1]].x] for item in idx], columns = ["ind1","y1","x1","ind2","y2","x2"])

        def check_expansion(df):
            """
                Drop overturning events that don't have a longitudinal expansion larger than min_exp
            """
            #calculate the longitudinal expansion
            check_exp = [(max(df_contour[int(row.ind1):int(row.ind2)+1].x) - min(df_contour[int(row.ind1):int(row.ind2)+1].x)) > min_exp for index,row in df.iterrows()]
            return df[check_exp]

        #apply all four routines and stop if no overturning events are left after a routine
        routines = [check_duplicates, check_joining, check_distance, check_expansion]
        i = -1
        while len(df_ot.index) > 0:
            i+=1
            df_ot = routines[i](df_ot).reset_index(drop=True)
            if i == 3:
                break


        #return the result in a DataFrame  
        if df_ot.empty:
            return pd.DataFrame([])
        else:
            def ot_to_grid(df):
                """
                     Get all grid cells that are enclosed by the rectangle describing an overturning event
                """
                #map the overturning events on the original grid, represented by a rectangle around the event
                temp = df_contour[int(df.ind1):int(df.ind2)+1]
                x, y = np.meshgrid(range(int(min(temp.x)), int(max(temp.x))+1),range(int(min(temp.y)), int(max(temp.y))+1)) 
                x, y = x.flatten(), y.flatten()
                return np.vstack((y,x)).T 
            
            overturnings_grids = [pd.DataFrame(ot_to_grid(row), columns = ["y", "x"])  for index,row in df_ot.iterrows()]
            return self.get_events(contour.date, overturnings_grids)
        
    def get_events(self, date, lodf):
        """
            This is an internal function to calculate several properties of the identified events.
        """
        def event_properties(df):
            """
                Calculate properties for each event.
            """
            #select time step and expand field for periodicity
            ds = self.dataset.sel({self._time_name:date})
            ds = xr.concat([ds, ds.isel({self._longitude_name:slice(0,self._periodic_add)})], dim=self._longitude_name)

            #select grid cells of a specific event and calculate several properties
            temp = ds.isel({self._latitude_name:df.y.to_xarray(), self._longitude_name:df.x.to_xarray()})
            area = np.round(np.sum(temp.weight.values),2)
            mean_var = np.round(np.average(temp[self._contour_variable].values),2)
            com = np.round(np.sum(df.multiply((temp.weight*temp[self._contour_variable]).values, axis = 0), axis = 0)/sum((temp.weight*temp[self._contour_variable]).values),2)
            
            #map coordinates of the event to the original grid
            df.x = self.dataset[self._longitude_name].values[df.x.values.astype("int")%self.dims[self._longitude_name]]
            df.y = self.dataset[self._latitude_name].values[df.y.values.astype("int")%self.dims[self._latitude_name]]
            df_coords = list(map(tuple,df.values))

            com[self._longitude_name]  = self.dataset[self._longitude_name].values[(np.round(com.x)%self.dims[self._longitude_name]).astype("int")]
            com[self._latitude_name]  = self.dataset[self._latitude_name].values[(np.round(com.y)%self.dims[self._latitude_name]).astype("int")]

            #if "mflux" is provided, calculate the intensity
            if "mflux" in self.variables:
                intensity = np.round(np.sum((temp.weight*temp["mflux"]).values)/area,2)
                return [df_coords, mean_var, area, intensity, (com[self._latitude_name] ,com[self._longitude_name])]
            else:
                return [df_coords, mean_var, area, (com[self._latitude_name] ,com[self._longitude_name])]

        #store properties in a DataFrame
        properties = [event_properties(item) for item in lodf]

        if "mflux" in self.variables:
            return pd.DataFrame([{"date":date, "coordinates":item[0], "mean_var":item[1], "area":item[2], "intensity":item[3], "com":item[4]} for item in properties])
        else:
            return pd.DataFrame([{"date":date, "coordinates":item[0], "mean_var":item[1], "area":item[2], "com":item[3]} for item in properties])
        
        
# ============================================================================
# post-processing functions
# ============================================================================


    def to_xarray(self, events, name = "flag"):
        
        """
        Create xarray.Dataset from events stored in a xarray.DataFrame
        Grid cells where an event is present are flagged with the value 1
        
        Parameters
        ----------
            events : pd.DataFrame
                DataFrame with the date and coordinates for every identified event
            name : string
                name of the xarray variable that is created
                
        Returns
        -------
            xarray.Dataset: int
                Dataset where the coordinates of the events are flagged with the value 1
        
        """
        
        logger.info('Creating xarray variable...')
        
        if events.empty:
            errmsg = 'Your events DataFrame is empty!'
            raise ValueError(errmsg)
        
        #get coordiantes of all events at the same time step
        events_concat = events.groupby(events.date).apply(lambda s: pd.concat([pd.DataFrame(row.coordinates, columns = ["lat", "lon"]) for index,row in s.iterrows()])).reset_index().drop("level_1", axis =1)
        
        #create emtpy Dataset with the same dimension as the original Dataset
        self.dataset[name] = xr.zeros_like(xr.DataArray(np.zeros((self.dims[self._time_name],self.dims[self._latitude_name],self.dims[self._longitude_name])), dims = [self._time_name, self._latitude_name, self._longitude_name]))
        
        #set coordinates of the events to the value 1
        self.dataset[name].loc[{self._time_name:events_concat.date.to_xarray(), self._latitude_name:events_concat[self._latitude_name].to_xarray(), self._longitude_name:events_concat[self._longitude_name].to_xarray()}] = np.ones(len(events_concat))
        self.dataset[name] = self.dataset[name].astype("int8")
        
        logger.info("Variable created: {}".format(name))
        
    def event_tracking(self, events, box = False):
        
        """
        Temporal tracking of events. 
        Events receive the same label if an event at step t spacially overlaps with an event at step t-1.
        
        Parameters
        ----------
            events : pd.DataFrame
                DataFrame with the date and coordinates for every identified event
            box : bool, optional
                If True, a rectangular box is used for the temporal tracking
                
        Returns
        -------
            xarray.Dataset: int
                Dataset where the coordinates of the events are flagged with the value 1
        
        """
        
        if events.empty:
            errmsg = 'Your events DataFrame is empty!'
            raise ValueError(errmsg)
            
        logger.info('Tracking events over time...')
        
        #step up progress bar
        dates = [pd.Timestamp(date).strftime('%Y-%m-%dT%H') for date in self.dataset[self._time_name].values]
        progress = tqdm(range(0,len(dates)-1), leave = True, position = 0)

        combine = []
        if box == True:
            #calculate rectengular boxes for each event
            boxes = [pd.DataFrame([(y,x) for y,x in itertools.product(set(np.asarray(row.coordinates)[:,0]), set(np.asarray(row.coordinates)[:,1]))], columns = ["lat", "lon"]) for index, row in events.iterrows()]
            for i,r in zip(range(0,len(dates)-1), progress):
                #check if events at step t overlap with events of t-1
                index_combinations = itertools.combinations(events[events.date.isin(dates[i:i+2])].index, r=2)
                combine.append([combination for combination in index_combinations if not set(map(tuple,boxes[events.index.get_loc(combination[0])].values)).isdisjoint(set(map(tuple,boxes[events.index.get_loc(combination[1])].values)))])
        else:
            #get coordinates of each event
            coords = [pd.DataFrame(row.coordinates, columns = ["lat", "lon"]) for index, row in events.iterrows()]
            for i,r in zip(range(0,len(dates)-1), progress):
                #check if events at step t overlap with events of t-1
                index_combinations = itertools.combinations(events[events.date.isin(dates[i:i+2])].index, r=2)
                combine.append([combination for combination in index_combinations if not set(map(tuple,coords[events.index.get_loc(combination[0])].values)).isdisjoint(set(map(tuple,coords[events.index.get_loc(combination[1])].values)))])

        #combine tracked indices to groups
        combine = list(itertools.chain.from_iterable(combine))
        combine = self.combine_shared(combine)

        events["label"] = events.index

        #assign lables to the events
        for item in combine:
            events.loc[item, "label"]= min(item)
            
        #select smallest possible index number for all events
        label = events.label.copy()
        for i in np.arange(len(set(events.label))):
            label[events.label == sorted(set(events.label))[i]] = i
        events.label = label
            
        #sort list by label and date
        self.labeled_events = events.sort_values(by=["label", "date"])
        
        
# ============================================================================
# plotting functions
# ============================================================================
        
    
    def plot_clim(self, variable, seasons = None, proj = "PlateCarree", smooth_passes = 10, periodic = True, labels = True, levels = None, cmap = None, title = ""):
        
        """
        Creates a simple climatological plot showing the occurrence frequency of the detected events. 
        
        Parameters
        ----------
            variable: string
                name of the xarray.Dataset variable containing the locations of the events flagged with the value 1
            seasons: list or array, optional
                months of the seasons of which the occurrence frequency should be calculated (e.g. [11,12,1])
            proj: string, optional
                name of cartopy projection
            smooth_passes: int or float, optional
                number of smoothing passes of the 5-point smoothing of the occurrence frequencies
            periodic: bool, optional
                If True, the first longitude is added at the end to close the gap in a polar projection
            labels: bool, optional
                If False, no labels are added to the plot
            levels: list or array, optional
                Colorbar levles. If not provided, default levels are used.
            cmap: string, optional
                Name of a valid cmap. If not provided, a default cmap is used. 
            title: string, optional
                Title of the plot
                
        Returns
        -------
            matplotlib.pyplot:
                Climatological plot of the occurrence frequencies.
        
        """
        
        if variable not in self.variables:
            errmsg = 'Variable {} not available! Use to_xarray() to transform the events to an xarray variable!'.format(variable)
            raise ValueError(errmsg)
        
        #define cartopy projections
        data_crs = ccrs.PlateCarree()
        proj = getattr(ccrs, proj)()

        #initialize figure
        fig, ax = plt.subplots(1,1, subplot_kw=dict(projection=proj), figsize=(12,8))

        #calculate occurrence frequencies, if provided for seasons
        if seasons is None:
            freq = xr.where(self.dataset[variable]>0,1,0).sum(dim=self._time_name)/self.dims[self._time_name]*100
        else:
            ds_season = self.dataset[variable].sel({self._time_name:self.dataset[self._time_name].dt.month.isin(seasons)})
            freq = xr.where(ds_season>0,1,0).sum(dim=self._time_name)/len(ds_season[self._time_name])*100

        #smooth occurrence frequencies
        if smooth_passes > 0:
            temp = freq.pad({self._longitude_name:5}, mode="wrap")
            temp = wrf.smooth2d(temp,smooth_passes)
            freq = temp.isel({self._longitude_name:slice(5,-5)})
        
        #add longitude to ensure that is no gap in a periodic field
        if periodic == True:
            freq = xr.concat([freq, freq.isel({self._longitude_name:0}).assign_coords({"lon":freq[self._longitude_name].max()+1})], dim = self._longitude_name)
            
        #define levels
        if levels is None:
            max_level = np.round(freq.max()*2)/2
            levels = set(np.arange(0, max_level, np.round(max_level/7*2)/2).tolist())
        
        #define cmap
        if cmap is None:
            def truncate_colormap(cmap, minval=0, maxval=1.0, n=100):
                new_cmap = colors.LinearSegmentedColormap.from_list(
                    'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                    cmap(np.linspace(minval, maxval, n)))
                return new_cmap
            cmap = truncate_colormap(plt.get_cmap("RdYlBu_r"), 0.3, 1)
            
        #plot frequencies
        p = freq.where(freq>0).plot.contourf(ax=ax, cmap = cmap, levels = levels,  transform=data_crs, add_colorbar = False, extend = "max")
        
        #define colorbar
        cax = fig.add_axes([ax.get_position().x1+0.05,ax.get_position().y0,0.015,ax.get_position().height])
        cbar = plt.colorbar(p, cax=cax, drawedges=True)
        cbar.ax.set_yticklabels(levels, fontsize=12, weight='bold')
        cbar.set_label(label="Occurrence frequency in %",size=12, fontweight="bold",labelpad=15)
        cbar.outline.set_color('black')
        cbar.outline.set_linewidth(2)
        cbar.dividers.set_color('black')
        cbar.dividers.set_linewidth(2)
        
        #add coast lines and grid lines
        ax.add_feature(cfeature.COASTLINE, color="dimgrey")
        gr = ax.gridlines(draw_labels=False, color="black", linestyle="dotted", linewidth = 1.1)
        
        #plot labels
        if labels == True:
            gr = ax.gridlines(draw_labels=True, color="black", linestyle="dotted", linewidth = 1.1)
            gr.xlabel_style = {'size': 12, 'color': 'black', "rotation":0}
            gr.ylabel_style = {'size': 12, 'color': 'black'}
        
        #make a circular cut out for the NorthPolarStereo projection
        if proj == ccrs.NorthPolarStereo():
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.add_patch(mpatches.Circle((0.5, 0.5), radius=0.5, color='k', linewidth=5, fill=False, transform = ax.transAxes))
            
        #set title
        ax.set_title(title, fontweight='bold',fontsize=20)
        
    
    def plot_step(self, flag_variable, variable, step, contour_level, proj = "PlateCarree",labels = True, levels = None, cmap = None, title = ""):
        
        """
        Creates a plot showing the events on one time step.
        
        Parameters
        ----------
            flag_variable: string
                name of the xarray.Dataset variable containing the locations of the events flagged with the value 1
            variable: string
                name of the xarray.Dataset variable that has been used to calculate the contours
            step: int or string
                index or name of a time step in the xarray.Dataset
            contour_level: list
                list of contour levels that are shown in the plot
            proj: string, optional
                name of cartopy projection
            labels: bool, optional
                If False, no labels are added to the plot
            levels: list or array, optional
                Colorbar levles. If not provided, default levels are used.
            cmap: string, optional
                Name of a valid cmap. If not provided, a default cmap is used. 
            title: string, optional
                Title of the plot
                
        Returns
        -------
            matplotlib.pyplot:
                Plot of one time step.
        
        """
        
        if variable not in self.variables:
            errmsg = 'Variable {} not available! Use to_xarray() to transform the events to an xarray variable!'.format(variable)
            raise ValueError(errmsg)
           
        #select date from data
        if (type(step) is str or type(step) is np.datetime64):
            ds = self.dataset.sel({self._time_name:step})
            ds = xr.concat([ds, ds.isel({self._longitude_name:0}).assign_coords({"lon":ds[self._longitude_name].max()+1})], dim = self._longitude_name)
            date = pd.to_datetime(str(self.dataset.sel({self._time_name:step})[self._time_name].values)).strftime('%Y-%m-%d')  
        else:
            try: 
                ds = self.dataset.isel({self._time_name:step})
                ds = xr.concat([ds, ds.isel({self._longitude_name:0}).assign_coords({"lon":ds[self._longitude_name].max()+1})], dim = self._longitude_name)
                date = pd.to_datetime(str(self.dataset.isel({self._time_name:step})[self._time_name].values)).strftime('%Y-%m-%d')
            except:
                errmsg = '{} as input for step is not supported! Step must be either an index or an coordinate of your data!'.format(step)
                raise ValueError(errmsg)                      

        #define cartopy projections        
        data_crs = ccrs.PlateCarree()
        proj = getattr(ccrs, proj)()

        #initialize figure
        fig, ax = plt.subplots(1,1, subplot_kw=dict(projection=proj), figsize=(12,8))
        
        #define levels
        if levels is None:
            max_level = np.round(ds[variable].max()*2)/2
            levels = set(np.arange(0, max_level, np.round(max_level/7*2)/2).tolist())
        
        #define cmap
        if cmap is None:
            cmap = "Blues"
            
        #plot variable field, contour line and flag_variable
        p = ds[variable].plot.contourf(ax=ax, cmap = cmap, levels = levels,  transform=data_crs, add_colorbar = False, extend = "max", alpha = 0.8)
        c = ds[variable].plot.contour(ax=ax, levels = np.asarray(contour_level).tolist(),  transform=data_crs, linewidths = 2, colors = "black")
        f = ds[flag_variable].where(ds[flag_variable]>0).plot.contourf(ax=ax, colors = ["white", "gold"], levels = [0,0.5],  transform=data_crs, add_colorbar = False)

        #define colorbar
        cax = fig.add_axes([ax.get_position().x1+0.05,ax.get_position().y0,0.015,ax.get_position().height])
        cbar = plt.colorbar(p, cax=cax, drawedges=True)
        cbar.ax.set_yticklabels(levels, fontsize=11, weight='bold')
        
        if all(x in self.dataset[variable].attrs for x in ["units", "long_name"]):
            cbar.set_label(label = self.dataset[variable].long_name + " [" +self.dataset[variable].units + "]", size = 12, fontweight="bold",labelpad=15)
            
        cbar.outline.set_color('black')
        cbar.outline.set_linewidth(2)
        cbar.dividers.set_color('black')
        cbar.dividers.set_linewidth(2)
        
        #add coast lines and grid lines
        ax.add_feature(cfeature.COASTLINE, color="dimgrey")
        gr = ax.gridlines(draw_labels=False, color="black", linestyle="dotted", linewidth = 1.1)
        
        #plot labels
        if labels == True:
            gr = ax.gridlines(draw_labels=True, color="black", linestyle="dotted", linewidth = 1.1)
            gr.xlabel_style = {'size': 12, 'color': 'black', "rotation":0}
            gr.ylabel_style = {'size': 12, 'color': 'black'}
        
        #make a circular cut out for the NorthPolarStereo projection
        if proj == ccrs.NorthPolarStereo():
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.add_patch(mpatches.Circle((0.5, 0.5), radius=0.5, color='k', linewidth=5, fill=False, transform = ax.transAxes))
        
        #set title
        ax.set_title(title, fontweight='bold',fontsize=20)
        
        #add the date to the figure
        plt.text(0.99,0.99, "Date: " + str(date), fontsize = 12, fontweight = "bold", ha='right', va='top', transform=ax.transAxes)
        
        
    def plot_tracks(self, events, proj = "PlateCarree", labels = True, min_path = 0, plot_events = False, title = ""):
        
        """
        Creates a plot showing the tracks of the temporally tracked events. 
        
        Parameters
        ----------
            events: pd.DataFrame
                DataFrame with the date, coordinates and label of every identified event
            proj: string, optional
                name of cartopy projection
            labels: bool, optional
                If False, no labels are added to the plot
            min_path: int, optional
                Minimal number of time steps an event has to be tracked
            plot_events: bool, optional
                If True, the events are also plotted by a shaded area
            title: string, optional
                Title of the plot
                
        Returns
        -------
            matplotlib.pyplot: 
                Plot of the tracks
        
        """
    
        if "label" not in events.columns:
            errmsg = 'No tracked events detected. Plase track events first.'
            raise ValueError(errmsg)
        
        #define cartopy projections
        data_crs = ccrs.PlateCarree()
        proj = getattr(ccrs, proj)()

        #initialize figure
        fig, ax = plt.subplots(1,1, subplot_kw=dict(projection=proj), figsize=(12,8))
        
        #group event data by label and plot each path
        for name, group in events.groupby("label"):
            if len(group)>min_path:
        
                lons = np.asarray(group.com.tolist())[:,1]
                lats = np.asarray(group.com.tolist())[:,0]
                
                #plot start point of each path
                ax.scatter(lons[0],lats[0], s=14, zorder=10,facecolors='none', edgecolor='black', transform=data_crs)
                ax.plot(lons[0], lats[0], ".", color = "red", transform=data_crs, alpha =0.7)
                
                #plot the coordinates of the events
                if plot_events == True:
                    for index,row in group.iterrows():
                        ax.plot(np.asarray(row.coordinates)[:,1], np.asarray(row.coordinates)[:,0], ".", color = "black", transform = data_crs, alpha = 0.2, markersize = 3)
                
                max_lon = max(self.dataset.lon.values)+1
                min_lon = min(self.dataset.lon.values)
                
                #check if the path needs to be split due to a crossing of the date border
                diffs = abs(np.diff(lons))>self.dims["lon"]/2
                split = [[i-1,i] for i in np.where(diffs)[0]+1]
                no_split = [[i-1,i] for i in np.where(~diffs)[0]+1]
                
                #plot paths that do not need to be split
                for item in no_split:
                    ev_seg = np.asarray(group[item[0]:item[1]+1].com.tolist())
                    
                    ax.plot(ev_seg[:,1], ev_seg[:,0], "-", transform=data_crs, color = "black")
                    
                #plot paths that have to be split
                for item in split:
                    ev_seg = np.asarray(group[item[0]:item[1]+1].com.tolist())

                    #upstream split
                    if np.diff(ev_seg[:,1])<0:
                        lons_plot = [(ev_seg[0,1], max_lon), (min_lon, ev_seg[1,1])]
                        lon_diffs = np.diff(lons_plot)

                        m = np.diff(ev_seg[:,0])[0]/np.sum(lon_diffs)
                        lats_plot = [(ev_seg[0,0],ev_seg[0,0]+lon_diffs[0][0]*m),(ev_seg[1,0]-lon_diffs[1][0]*m, ev_seg[1,0])]

                        for lon, lat in zip(lons_plot, lats_plot):
                            ax.plot(lon, lat, "-", "-", transform=data_crs, color = "black")

                    #downstream split
                    else:
                        lons_plot = [(ev_seg[0,1], min_lon), (max_lon, ev_seg[1,1])]
                        lon_diffs = np.diff(lons_plot)

                        m = np.diff(ev_seg[:,0])[0]/np.sum(lon_diffs)
                        lats_plot = [(ev_seg[0,0],ev_seg[0,0]+lon_diffs[0][0]*m),(ev_seg[1,0]-lon_diffs[1][0]*m, ev_seg[1,0])]

                        for lon, lat in zip(lons_plot, lats_plot):
                            ax.plot(lon, lat, "-", "-", transform=data_crs, color = "black")        
                            
        #plot invisible weigths to get a plot of the full grid       
        p = self.dataset.weight.plot.contourf(ax=ax, transform=data_crs, add_colorbar = False, alpha =0)

        #add coast lines and grid lines
        ax.add_feature(cfeature.COASTLINE, color="dimgrey")
        gr = ax.gridlines(draw_labels=False, color="black", linestyle="dotted", linewidth = 1.1)

        #plot labels
        if labels == True:
            gr = ax.gridlines(draw_labels=True, color="black", linestyle="dotted", linewidth = 1.1)
            gr.xlabel_style = {'size': 12, 'color': 'black', "rotation":0}
            gr.ylabel_style = {'size': 12, 'color': 'black'}

        #make a circular cut out for the NorthPolarStereo projection
        if proj == ccrs.NorthPolarStereo():
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.add_patch(mpatches.Circle((0.5, 0.5), radius=0.5, color='k', linewidth=5, fill=False, transform = ax.transAxes))
        
        #set title
        ax.set_title(title, fontweight='bold',fontsize=20)
        
        
# ============================================================================
# utility functions
# ============================================================================
    
    
    def combine_shared(self,lst):
        
        """
        This is an internal function that combines all elements of a list that have at least one element in common.  
        
        Parameters
        ----------
            lst : list
                list of integer or float values
                
        Returns
        -------
            list: int or float
                list of combined elements
        
        """
        
        l = lst
        output = []
        while len(l)>0:
            first, *rest = l
            first = set(first)

            lf = -1
            while len(first)>lf:
                lf = len(first)

                rest2 = []
                for r in rest:
                    if len(first.intersection(set(r)))>0:
                        first |= set(r)
                    else:
                        rest2.append(r)     
                rest = rest2

            output.append(list(first))
            l = rest
            
        return output

# ============================================================================