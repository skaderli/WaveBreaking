#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""

Description
-----------
 
    add description...

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
from numpy.core import datetime64
from skimage import measure
from scipy import spatial
import itertools as itertools
from tqdm import tqdm
import wrf as wrf
from shapely.geometry import Point, LineString, MultiPoint, Polygon
import shapely.vectorized
from sklearn.metrics import DistanceMetric
dist = DistanceMetric.get_metric('haversine')

# plotting
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import seaborn as sns
import matplotlib.colors as colors

# logs
import logging
logger = logging.getLogger(__name__)
#logger.setLevel(logging.DEBUG)
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
            except IOError:
                raise IOError("Unkown fileformat. Known formats are netcdf or csv.")

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
            Empty contrack container.\n\
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

# ----------------------------------------------------------------------------
# Read / Import / Save data
    
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
            errmsg = 'contrack() is already set!'
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
 
# ----------------------------------------------------------------------------
# Set up / Check dimensions
   
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
            if isinstance(var, datetime64):
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

# ----------------------------------------------------------------------------
# field processing

    def calculate_momentum_flux(self, variable_zonal, variable_meridional):

        if hasattr(self, '_time_name'):
            logger.info("Check dimensions... OK")
            pass
        else:
            logger.info("Set-up dimensions... DONE")
            self.set_up()

        logger.info('Calculating momentum flux...')

        dzm_u = self.dataset[variable_zonal].resample({self._time_name:'1D'}).mean(self._longitude_name)
        dzm_v = self.dataset[variable_meridional].resample({self._time_name:'1D'}).mean(self._longitude_name)

        u_prime = self.dataset[variable_zonal]-dzm_u
        v_prime = self.dataset[variable_meridional]-dzm_v

        flux = u_prime*v_prime
        self.dataset["mflux"] = flux

        logger.info("Variable created: 'mflux'")
        
        
    def calculate_smoothed_field(self, variable, passes):
        
        if hasattr(self, '_time_name'):
            logger.info("Check dimensions... OK")
            pass
        else:
            logger.info("Set-up dimensions... DONE")
            self.set_up()

        logger.info('Calculating smoothed field...')

        temp = self.dataset[variable].pad({self._longitude_name:5}, mode="wrap")
        temp = wrf.smooth2d(temp,passes)

        self.dataset[temp.name] = temp.isel({self._longitude_name:slice(5,-5)})

        logger.info("Variable created: '{}'".format(temp.name))
    
# ----------------------------------------------------------------------------
# contour calculations
        
    def get_contours(self, variable, level, periodic_add = 120):
  
        if hasattr(self, '_time_name'):
            logger.info("Check dimensions... OK")
            pass
        else:
            logger.info("Set-up dimensions... DONE")
            self.set_up()
        
        logger.info('Calculating contour lines...')
        times = self.dataset[self._time_name]
        progress = tqdm(range(0,len(times)), leave = True, position = 0)
        
        contours = [self.get_contour_timestep(step, level=level, variable=variable, periodic_add=periodic_add) for step,p in zip(times,progress)]
        
        self.contours = pd.concat([item[0] for item in contours]).reset_index(drop=True)
        
        contours_wb_calc = pd.concat([item[1] for item in contours])
        self._contours_wb_calc = contours_wb_calc[contours_wb_calc.exp == contours_wb_calc.exp.max()].reset_index(drop=True)
        
        logger.info('{} contour(s) identified!'.format(len(self.contours)))
        
    def get_contour_timestep(self, step, variable, level, periodic_add):

        #define data
        ds = self.dataset.sel({self._time_name:step})[variable]
        #expand field for periodicity
        ds = xr.concat([ds, ds.isel({self._longitude_name:slice(0,periodic_add)})], dim=self._longitude_name)

        #get contours (indices in array coordinates)
        contours_from_measure = measure.find_contours(ds.values, level)

        #contours for wave breaking calculation
        contours_index_expanded = [pd.DataFrame(list(dict.fromkeys(map(tuple,np.round(item).astype("int")))), columns=["y", "x"]) for item in contours_from_measure]

        #original contours for output
        contours_index_original = [np.c_[item.y.values%self.dims[self._latitude_name], item.x.values%self.dims[self._longitude_name]] for item in contours_index_expanded]
        contours_index_original =  [list(dict.fromkeys(map(tuple,item))) for item in contours_index_original]
        
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
        
        contour_coords_original = [pd.DataFrame(np.c_[self.dataset[self._latitude_name].values[item[:,0].astype("int")], 
                                                      self.dataset[self._longitude_name].values[item[:,1].astype("int")]], 
                                                columns = [self._latitude_name, self._longitude_name]) for item in contours_index_original]
        
        #return contour coordinates in a DataFrame
        def contour_to_df(lodf):
            date = pd.to_datetime(str(step.values)).strftime('%Y-%m-%dT%H')
            x_exp = [item.iloc[:,1].max() - item.iloc[:,1].min() for item in lodf]
            mean_lat = [item.iloc[:,0].mean() for item in lodf]

            return pd.DataFrame([{"date":date, "coords":c, "level":level, "exp":x, "mean_lat":m} for c,x,m in zip(lodf, x_exp, mean_lat) if len(c)>1])

        return contour_to_df(contour_coords_original), contour_to_df(contours_index_expanded)

    
        
# ----------------------------------------------------------------------------
# index calculations        
        
    def get_streamers(self, variable, level, geo_dis = 800, cont_dis = 1500, periodic_add = 120):
        
        if hasattr(self, '_contours_wb_calc'):
            if self._contours_wb_calc.level.iloc[0] == level:
                logger.info("Check contours... OK")
                pass
            else:
                logger.info("Check contours... Available contours are on a different level! Calculating new contours...")
                self.get_contours(variable=variable, level=level, periodic_add = periodic_add)
        else:
            logger.info("Check contours... No contours found! Calculating contours...")
            self.get_contours(variable = variable, level=level, periodic_add = periodic_add)
            
        if "mflux" in self.variables:
            logger.info("Momentum flux variable detected. Intensity will be calculated.")
        else:
            logger.info("No momentum flux detected. Intensity will not be calculated.")
       
             
        # create variable "weight" representing the weight of each grid cell
        weight_lat = np.cos(self.dataset[self._latitude_name].values*np.pi/180)
        self.dataset["weight"] = xr.DataArray(np.ones((self.dims[self._latitude_name], self.dims[self._longitude_name])) * np.array((111 * 1 * 111 * 1 * weight_lat)).astype(np.float32)[:, None], dims = [self._latitude_name, self._longitude_name])
        
        logger.info('Calculating streamers...')
        times = self.dataset[self._time_name]
        progress = tqdm(range(0,len(times)), leave = True, position = 0)
        
        self.streamers = pd.concat([self.get_streamers_timestep(row, variable, geo_dis, cont_dis, periodic_add) for (index,row),p in zip(self._contours_wb_calc.iterrows(),progress)]).reset_index(drop=True)
        
        logger.info('{} streamer(s) identified!'.format(len(self.streamers)))
        
    def get_streamers_timestep(self, contour, variable, geo_dis, cont_dis, periodic_add):
        
        #In a first step, calculate all possible basepoints    
        df_contour = contour.coords
        contour_line = LineString([Point(row.x, row.y) for index, row in df_contour.iterrows()]) 

        #calculate geographic distance between all contour coordinates
        geo_matrix = np.triu((dist.pairwise((np.radians(df_contour.loc[:,["y", "x"]])))*6371))

        #calculate the connecting distance between all contour coordinates
        on = np.insert(np.diagonal(geo_matrix,1),0,0)
        on_mat = np.triu(np.tile(on,(len(on),1)), k=1)
        cont_matrix = np.cumsum(on_mat, axis = 1)

        #get indices of all contour coordinates that fulfill both conditions
        check = (geo_matrix<geo_dis) * (cont_matrix>cont_dis)
        check = np.transpose(check.nonzero())

        #store all coordinates in a DataFrame
        df1 = df_contour[["y","x"]].iloc[check[:,0]].reset_index().rename(columns = {"y":"y1", "x":"x1", "index":"ind1"})
        df2 = df_contour[["y","x"]].iloc[check[:,1]].reset_index().rename(columns = {"y":"y2", "x":"x2", "index":"ind2"})
        df_bp = pd.concat([df1, df2], axis = 1)
        df_bp = df_bp.drop(df_bp[np.abs(df_bp.x1-df_bp.x2)>120].index)

        #In second step, several routines are applied to reduce the number of basepoints
        def check_duplicates(df):
            temp = pd.concat([df.y1, df.x1%self.dims[self._longitude_name], df.y2, df.x2%self.dims[self._longitude_name]], axis = 1)
            temp = temp[temp.duplicated(keep=False)]
            if len(temp) == 0:
                check = []
            else:
                check = list(np.asarray(temp.groupby(list(temp)).apply(lambda x: tuple(x.index)).tolist())[:,1])

            return df.drop(check)

        def check_intersections(df):
            #check for intersections with the contour line 
            df["dline"] = [LineString([(row.x1, row.y1),(row.x2, row.y2)]) for index, row in df.iterrows()]
            check_intersections = [row.dline.touches(contour_line) for index, row in df.iterrows()]

            return df[check_intersections]

        def check_overlapping(df):
            #check for overlapping of the contour segment described by basepoints
            ranges = [range(int(r2.ind1),int(r2.ind2+1)) for i2,r2 in df.iterrows()]
            index_combinations = list(itertools.permutations(df.index, r=2))
            check_overlapping = set([item[0] for item in index_combinations if (ranges[item[0]][0] in ranges[item[1]] and ranges[item[0]][-1] in ranges[item[1]])])

            return df.drop(check_overlapping)

        def check_groups(df):
            index_combinations = np.asarray(list(itertools.combinations_with_replacement(df.index, r=2)))
            check_crossing = [df.iloc[item[0]].dline.intersects(df.iloc[item[1]].dline) for item in index_combinations]

            groups = self.combine_shared(index_combinations[check_crossing])
            keep_index = [item[np.argmax([on[int(row.ind1):int(row.ind2)+1].sum() for index, row in df.iloc[item].iterrows() ])] for item in groups]

            return df.iloc[keep_index]

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
                x, y = np.meshgrid(np.arange(0,self.dims[self._longitude_name]+periodic_add), np.arange(0,self.dims[self._latitude_name]))
                x, y = x.flatten(), y.flatten()
                mask = shapely.vectorized.contains(Polygon(df_contour[int(df.ind1):int(df.ind2)+1]),y,x)
                return np.r_[df_contour[int(df.ind1):int(df.ind2)+1].values, np.c_[y,x][mask]]


        streamer_grids = [pd.DataFrame(bp_to_grid(row), columns = ["y", "x"])  for index,row in df_bp.iterrows()]
        return self.get_events(contour.date, streamer_grids, variable, periodic_add)
    
    def get_overturnings(self, variable, level, periodic_add = 120):
        
        if hasattr(self, '_contours_wb_calc'):
            if self._contours_wb_calc.level.iloc[0] == level:
                logger.info("Check contours... OK")
                pass
            else:
                logger.info("Check contours... Available contours are on a different level! Calculating new contours...")
                self.get_contours(variable=variable, level=level, periodic_add = periodic_add)
        else:
            logger.info("Check contours... No contours found! Calculating contours...")
            self.get_contours(variable = variable, level=level, periodic_add = periodic_add)
            
        if "mflux" in self.variables:
            logger.info("Momentum flux variable detected. Intensity will be calculated.")
        else:
            logger.info("No momentum flux detected. Intensity will not be calculated.")
       
             
        # create variable "weight" representing the weight of each grid cell
        weight_lat = np.cos(self.dataset[self._latitude_name].values*np.pi/180)
        self.dataset["weight"] = xr.DataArray(np.ones((self.dims[self._latitude_name], self.dims[self._longitude_name])) * np.array((111 * 1 * 111 * 1 * weight_lat)).astype(np.float32)[:, None], dims = [self._latitude_name, self._longitude_name])
        
        logger.info('Calculating overturning events...')
        times = self.dataset[self._time_name]
        progress = tqdm(range(0,len(times)), leave = True, position = 0)
        
        self.overturnings = pd.concat([self.get_overturnings_timestep(row, variable, periodic_add) for (index,row),p in zip(self._contours_wb_calc.iterrows(),progress)]).reset_index(drop=True)
        
        logger.info('{} overturning event(s) identified!'.format(len(self.overturnings)))
    
    def get_overturnings_timestep(self, contour, variable, periodic_add):
        
        #In a first step, calculate all overturning longitudes 
        df_contour = contour.coords   
        lons, counts = np.unique(df_contour.x, return_counts = True)

        ot_lons = lons[counts >= 3]
        ot_borders = [(df_contour[df_contour.x == lon].index[0], df_contour[df_contour.x == lon].index[-1]) for lon in ot_lons]

        df_ot =  pd.DataFrame([[item[0], df_contour.iloc[item[0]].y, df_contour.iloc[item[0]].x, 
                            item[1],  df_contour.iloc[item[1]].y, df_contour.iloc[item[1]].x] for item in ot_borders], columns = ["ind1","y1","x1","ind2","y2","x2"])

        #In second step, several routines are applied to reduce the number of events

        def check_duplicates(df):
            temp = pd.concat([df.y1, df.x1%self.dims[self._longitude_name], df.y2, df.x2%self.dims[self._longitude_name]], axis = 1)
            temp = temp[temp.duplicated(keep=False)]
            if len(temp) == 0:
                check = []
            else:
                check = list(np.asarray(temp.groupby(list(temp)).apply(lambda x: tuple(x.index)).tolist())[:,1])

            return df.drop(check)

        def check_joining(df):
            #check for overlapping of the contour segment described by basepoints
            ranges = [range(int(r2.ind1),int(r2.ind2)+1) for i2,r2 in df.iterrows()]
            index_combinations = np.asarray(list(itertools.combinations_with_replacement(df.index, r=2)))
            check_join = [not elem for elem in [set(ranges[item[0]]).isdisjoint(ranges[item[1]]) for item in index_combinations]]

            groups = self.combine_shared(index_combinations[check_join])

            idx = [(df_ot.iloc[item].ind1.min(), df_ot.iloc[item].ind2.max()) for item in groups]

            ret = pd.DataFrame([[item[0], df_contour.iloc[item[0]].y, df_contour.iloc[item[0]].x, 
                                 item[1], df_contour.iloc[item[1]].y, df_contour.iloc[item[1]].x] for item in idx], columns = ["ind1","y1","x1","ind2","y2","x2"])
            return ret

        def check_distance(df):
            geo_matrix = np.triu((dist.pairwise(np.radians(df.loc[:,['y2','x2']]), np.radians(df.loc[:,['y1','x1']]))*6371))
            indices_check_dis = [(i,j) for i,j in itertools.combinations_with_replacement(df.index, r=2) if geo_matrix[i,j] < 500]

            groups = self.combine_shared(indices_check_dis+[(i,i) for i in df.index])

            idx = [(df.iloc[item].ind1.min(), df.iloc[item].ind2.max()) for item in groups]

            ret = pd.DataFrame([[item[0], df_contour.iloc[item[0]].y, df_contour.iloc[item[0]].x, 
                                 item[1], df_contour.iloc[item[1]].y, df_contour.iloc[item[1]].x] for item in idx], columns = ["ind1","y1","x1","ind2","y2","x2"])
            return ret

        def check_expansion(df):
            check_exp = [(max(df_contour[int(row.ind1):int(row.ind2)+1].x) - min(df_contour[int(row.ind1):int(row.ind2)+1].x)) > 5 for index,row in df.iterrows()]
            return df[check_exp]

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
                temp = df_contour[int(df.ind1):int(df.ind2)+1]
                x, y = np.meshgrid(range(int(min(temp.x)), int(max(temp.x))+1),range(int(min(temp.y)), int(max(temp.y))+1)) 
                x, y = x.flatten(), y.flatten()
                return np.vstack((y,x)).T 
            
        overturnings_grids = [pd.DataFrame(ot_to_grid(row), columns = ["y", "x"])  for index,row in df_ot.iterrows()]
        return self.get_events(contour.date, overturnings_grids, variable, periodic_add)
        
    def get_events(self, date, lodf, variable, periodic_add):
        
        def event_properties(df):
            
            ds = self.dataset.sel({self._time_name:date})
            ds = xr.concat([ds, ds.isel({self._longitude_name:slice(0,periodic_add)})], dim=self._longitude_name)

            temp = ds.isel({self._latitude_name:df.y.to_xarray(), self._longitude_name:df.x.to_xarray()})
            area = np.round(np.sum(temp.weight.values),2)
            mean_var = np.round(np.average(temp[variable].values),2)
            com = np.round(np.sum(df.multiply((temp.weight*temp[variable]).values, axis = 0), axis = 0)/sum((temp.weight*temp[variable]).values),2)
            
            df.x = self.dataset[self._longitude_name].values[df.x.values.astype("int")%self.dims[self._longitude_name]]
            df.y = self.dataset[self._latitude_name].values[df.y.values.astype("int")%self.dims[self._latitude_name]]
            df = df.rename(columns = {"x":self._longitude_name, "y":self._latitude_name})

            com[self._longitude_name]  = self.dataset[self._longitude_name].values[(np.round(com.x)%self.dims[self._longitude_name]).astype("int")]
            com[self._latitude_name]  = self.dataset[self._latitude_name].values[(np.round(com.y)%self.dims[self._latitude_name]).astype("int")]

            if "mflux" in self.variables:
                intensity = np.round(np.sum((temp.weight*temp["mflux"]).values)/area,2)
                return [df, mean_var, area, intensity, (com[self._latitude_name] ,com[self._longitude_name])]
            else:
                return [df, mean_var, area, (com[self._latitude_name] ,com[self._longitude_name])]

        properties = [event_properties(item) for item in lodf]

        if "mflux" in self.variables:
            return pd.DataFrame([{"date":date, "coords":item[0], "mean_var":item[1], "area":item[2], "intensity":item[3], "com":item[4]} for item in properties])
        else:
            return pd.DataFrame([{"date":date, "coords":item[0], "mean_var":item[1], "area":item[2], "com":item[3]} for item in properties])
        
# ----------------------------------------------------------------------------
# post-processing and plotting functions

    def to_xarray(self, events, name = "flag"):
        
        logger.info('Creating xarray variable...')
        
        if events.empty:
            errmsg = 'Your events DataFrame is empty!'
            raise ValueError(errmsg)
        
        events_concat = events.groupby(events.date).apply(lambda s: pd.concat([row.coords for index,row in s.iterrows()])).reset_index().drop("level_1", axis =1)
        self.dataset[name] = xr.zeros_like(xr.DataArray(np.zeros((self.dims[self._time_name],self.dims[self._latitude_name],self.dims[self._longitude_name])), dims = [self._time_name, self._latitude_name, self._longitude_name]))
        self.dataset[name].loc[{self._time_name:events_concat.date.to_xarray(), self._latitude_name:events_concat[self._latitude_name].to_xarray(), self._longitude_name:events_concat[self._longitude_name].to_xarray()}] = np.ones(len(events_concat))
        self.dataset[name] = self.dataset[name].astype("int8")
        
        logger.info("Variable created: {}".format(name))
        
    def event_tracking(self, events):
        
        if events.empty:
            errmsg = 'Your events DataFrame is empty!'
            raise ValueError(errmsg)
            
        logger.info('Tracking events over time...')
        dates = [pd.Timestamp(date).strftime('%Y-%m-%d') for date in self.dataset[self._time_name].values]
        progress = tqdm(range(0,len(dates)-1), leave = True, position = 0)

        combine = []

        for i,r in zip(range(0,len(dates)-1), progress):
            index_combinations = itertools.combinations(events[events.date.isin(dates[i:i+2])].index, r=2)
            combine.append([combination for combination in index_combinations if not set(map(tuple,events.loc[combination[0]].coords.values)).isdisjoint(set(map(tuple,events.loc[combination[1]].coords.values)))])

        combine = list(itertools.chain.from_iterable(combine))
        combine = self.combine_shared(combine)

        events["label"] = events.index

        for item in combine:
            events.loc[item, "label"]= min(item)
            
        label = events.label.copy()
        for i in np.arange(len(set(events.label))):
            label[events.label == sorted(set(events.label))[i]] = i
            
        events.label = label
            
        self.labeled_events = events.sort_values(by=["label"])
        
    def plot_clim(self, variable, seasons = None, proj = ccrs.PlateCarree(), smooth_passes = 10, periodic = True, labels = True, levels = None, cmap = None, title = ""):
        
        if variable not in self.variables:
            errmsg = 'Variable {} not available! Use to_xarray() to transform the events to an xarray variable!'.format(variable)
            raise ValueError(errmsg)
        
        data_crs = ccrs.PlateCarree()
        proj = proj

        fig, ax = plt.subplots(1,1, subplot_kw=dict(projection=proj), figsize=(12,8))

        if seasons is None:
            freq = xr.where(self.dataset[variable]>0,1,0).sum(dim=self._time_name)/self.dims[self._time_name]*100
        else:
            ds_season = self.dataset[variable].sel({self._time_name:self.dataset[self._time_name].dt.month.isin(seasons)})
            freq = xr.where(ds_season>0,1,0).sum(dim=self._time_name)/len(ds_season[self._time_name])*100

        if smooth_passes > 0:
            temp = freq.pad({self._longitude_name:5}, mode="wrap")
            temp = wrf.smooth2d(temp,smooth_passes)
            freq = temp.isel({self._longitude_name:slice(5,-5)})
        
        if periodic == True:
            freq = xr.concat([freq, freq.isel({self._longitude_name:0}).assign_coords({"lon":freq[self._longitude_name].max()+1})], dim = self._longitude_name)
            
        if levels is None:
            max_level = np.round(freq.max()*2)/2
            levels = set(np.arange(0, max_level, np.round(max_level/7*2)/2).tolist())
        
        # define cmap
        if cmap is None:
            def truncate_colormap(cmap, minval=0, maxval=1.0, n=100):
                new_cmap = colors.LinearSegmentedColormap.from_list(
                    'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                    cmap(np.linspace(minval, maxval, n)))
                return new_cmap
            cmap = truncate_colormap(plt.get_cmap("RdYlBu_r"), 0.3, 1)
            
        p = freq.where(freq>0).plot.contourf(ax=ax, cmap = cmap, levels = levels,  transform=data_crs, add_colorbar = False, extend = "max")
        
        cax = fig.add_axes([ax.get_position().x1+0.05,ax.get_position().y0,0.015,ax.get_position().height])
        cbar = plt.colorbar(p, cax=cax, drawedges=True)
        cbar.ax.set_yticklabels(levels, fontsize=12, weight='bold')
        cbar.set_label(label="Occurrence frequency in %",size=12, fontweight="bold",labelpad=15)

        cbar.outline.set_color('black')
        cbar.outline.set_linewidth(2)

        cbar.dividers.set_color('black')
        cbar.dividers.set_linewidth(2)
        
        ax.add_feature(cfeature.COASTLINE, color="dimgrey")
        gr = ax.gridlines(draw_labels=False, color="black", linestyle="dotted", linewidth = 1.1)
        
        if labels == True:
            gr = ax.gridlines(draw_labels=True, color="black", linestyle="dotted", linewidth = 1.1)
            gr.xlabel_style = {'size': 12, 'color': 'black', "rotation":0}
            gr.ylabel_style = {'size': 12, 'color': 'black'}
        
        if proj == ccrs.NorthPolarStereo():
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.add_patch(mpatches.Circle((0.5, 0.5), radius=0.5, color='k', linewidth=5, fill=False, transform = ax.transAxes))
            
        ax.set_title(title, fontweight='bold',fontsize=20)
        
    
    def plot_step(self, flag_variable, variable, step, contour_level, proj = ccrs.PlateCarree(), levels = None, labels = True, cmap = None, title = ""):
        
        if variable not in self.variables:
            errmsg = 'Variable {} not available! Use to_xarray() to transform the events to an xarray variable!'.format(variable)
            raise ValueError(errmsg)
                        
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

        data_crs = ccrs.PlateCarree()
        proj = proj

        fig, ax = plt.subplots(1,1, subplot_kw=dict(projection=proj), figsize=(12,8))
        
        if levels is None:
            max_level = np.round(ds[variable].max()*2)/2
            levels = set(np.arange(0, max_level, np.round(max_level/7*2)/2).tolist())
        
        if cmap is None:
            cmap = "Blues"
            
        p = ds[variable].plot.contourf(ax=ax, cmap = cmap, levels = levels,  transform=data_crs, add_colorbar = False, extend = "max", alpha = 0.8)
        c = ds[variable].plot.contour(ax=ax, levels = np.asarray(contour_level).tolist(),  transform=data_crs, linewidths = 2, colors = "black")
        f = ds[flag_variable].where(ds[flag_variable]>0).plot.contourf(ax=ax, colors = ["white", "gold"], levels = [0,0.5],  transform=data_crs, add_colorbar = False)

        cax = fig.add_axes([ax.get_position().x1+0.05,ax.get_position().y0,0.015,ax.get_position().height])
        cbar = plt.colorbar(p, cax=cax, drawedges=True)
        cbar.ax.set_yticklabels(levels, fontsize=11, weight='bold')
        cbar.set_label(label = self.dataset[variable].long_name + " [" +self.dataset[variable].units + "]", size = 12, fontweight="bold",labelpad=15)

        cbar.outline.set_color('black')
        cbar.outline.set_linewidth(2)

        cbar.dividers.set_color('black')
        cbar.dividers.set_linewidth(2)

        ax.add_feature(cfeature.COASTLINE, color="dimgrey")
        gr = ax.gridlines(draw_labels=False, color="black", linestyle="dotted", linewidth = 1.1)
        
        if labels == True:
            gr = ax.gridlines(draw_labels=True, color="black", linestyle="dotted", linewidth = 1.1)
            gr.xlabel_style = {'size': 12, 'color': 'black', "rotation":0}
            gr.ylabel_style = {'size': 12, 'color': 'black'}
        
        if proj == ccrs.NorthPolarStereo():
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.add_patch(mpatches.Circle((0.5, 0.5), radius=0.5, color='k', linewidth=5, fill=False, transform = ax.transAxes))
            
        ax.set_title(title, fontweight='bold',fontsize=20)
        
        plt.text(0.99,0.99, "Date: " + str(date), fontsize = 12, fontweight = "bold", ha='right', va='top', transform=ax.transAxes)
        
        
    def plot_tracks(self, events, proj = ccrs.PlateCarree(), title = "", labels = True, min_path = 0):
    
        if "label" not in events.columns:
            errmsg = 'No tracked events detected. Plase track events first.'
            raise ValueError(errmsg)
            
        data_crs = ccrs.PlateCarree()
        proj = proj

        fig, ax = plt.subplots(1,1, subplot_kw=dict(projection=proj), figsize=(12,8))
        
        for name, group in events.groupby("label"):
            if len(group)>min_path:
                
                lons = np.asarray(group.com.tolist())[:,1]
                lats = np.asarray(group.com.tolist())[:,0]

                ax.plot(lons, lats, "-", transform=data_crs, color = "black")
                ax.scatter(lons[0],lats[0], s=14, zorder=10,facecolors='none', edgecolor='black', transform=data_crs)
                ax.plot(lons[0], lats[0], ".", color = "red", transform=data_crs, alpha =0.7)
                """
                else:
                    cut = np.where(np.diff(lons)<0)[0][0]+1
                    max_lon = max(wb.dataset.lon.values)
                    min_lon = min(wb.dataset.lon.values)
                    seg1 = max(wb.dataset.lon.values)-lons[cut-1]
                    seg2 = lons[cut]-min(wb.dataset.lon.values)
                    m = (lats[cut]-lats[cut-1])/(seg1+seg2+1)

                    lon_segs = [lons[:cut]+[max_lon+1], [min_lon]+lons[cut:]]
                    lat_segs = [lats[:cut]+[lats[cut-1]+(seg1+1)*m], [lats[cut-1]+(seg1+1)*m]+lats[cut:]]
                    for lon, lat in zip(lon_segs, lat_segs):
                        ax.plot(lon, lat, "-", transform=data_crs, color = "black")
                    ax.scatter(lon_segs[0][0],lat_segs[0][0], s=14, zorder=10,facecolors='none', edgecolor='black', transform=data_crs)
                    ax.plot(lon_segs[0][0], lat_segs[0][0], ".", color = "red", transform=data_crs, alpha =0.7)
                """
                    
        p = self.dataset.weight.plot.contourf(ax=ax, transform=data_crs, add_colorbar = False, alpha =0)

        ax.add_feature(cfeature.COASTLINE, color="dimgrey")
        gr = ax.gridlines(draw_labels=False, color="black", linestyle="dotted", linewidth = 1.1)

        if labels == True:
            gr = ax.gridlines(draw_labels=True, color="black", linestyle="dotted", linewidth = 1.1)
            gr.xlabel_style = {'size': 12, 'color': 'black', "rotation":0}
            gr.ylabel_style = {'size': 12, 'color': 'black'}

        if proj == ccrs.NorthPolarStereo():
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.add_patch(mpatches.Circle((0.5, 0.5), radius=0.5, color='k', linewidth=5, fill=False, transform = ax.transAxes))
            
        ax.set_title(title, fontweight='bold',fontsize=20)
        
        
        
        
# ----------------------------------------------------------------------------
# utility functions
    
    def combine_shared(self,lst):
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
