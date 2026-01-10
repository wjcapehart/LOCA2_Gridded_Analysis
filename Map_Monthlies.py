#!/usr/bin/env python
# coding: utf-8

# # Plot Maps

# In[1]:


##########################################################
#
# Library Calls.
#

# loading numpy

import numpy             as np

# loading matplotlib

import matplotlib.pyplot as plt

# loading xarray

import xarray            as xr
import h5netcdf          as h5netcdf

import datetime as datetime

# Loading pandas

import pandas            as pd

import os as os

import subprocess as subprocess


import matplotlib.pyplot as plt
import matplotlib        as mpl


pd.set_option('display.max_rows', None)

#
##########################################################


# In[ ]:





# In[2]:


##########################################################
#
# Time Coordinates
#

years_start = np.arange(start =    1950, 
                        stop  =    2071,   dtype=np.int32)
years_end   = np.array(years_start + 29, dtype=np.int32) 

years_middle = np.array(years_start/2 + years_end/2, dtype=np.float32)
years_bounds = np.array([years_start,years_end]).transpose()
n_runnings = len(years_start)

climatology_bounds = np.array([years_start,years_end]).transpose()



my_years = [1990.5, 2050.5, 2084.5]

yearall = xr.DataArray(name = "year",
                      data   = years_middle,
                      coords = {"year": years_middle},
                      attrs  = {"description": "middle calendar year for 30-yr period",
                                "long_name": "middle calendar year for 30-yr period",
                                "climatology_bounds":"climatology_bounds"})


year_bnds = xr.DataArray(name = "climatology_bounds",
                      data   =  years_bounds,
                      coords = {"year": yearall,
                                "bnds": 2},
                      attrs  = {"description": "boundary calendar years for 30-yr period",
                                "long_name": "boundary calendar years for 30-yr period"})


year_bnds

# year_bnds.to_netcdf("./climatological_bounds.nc")
#
##########################################################


# In[ ]:





# In[3]:


scenarios            = np.array(["historical", 
                                 "ssp245", 
                                 "ssp370", 
                                 "ssp585"])

scenario = xr.DataArray(name   = "scenario",
                        coords = {"scenario":scenarios},
                        data   = scenarios,
                        attrs  = {  "long_name":"scenario",
                                  "description":"scenario"})
scenario


# In[ ]:





# In[4]:


##########################################################
#
# Files
#

root_data_dir = "/data/DATASETS/LOCA_MACA_Ensembles/LOCA2/LOCA2_CONUS/Climate_CONUS/Monthly/"


mf_tasmax = [root_data_dir + 'historical/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmax___ALLRANK01___historical.nc',
             root_data_dir +     'ssp245/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmax___ALLRANK01___ssp245.nc',
             root_data_dir +     'ssp370/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmax___ALLRANK01___ssp370.nc',
             root_data_dir +     'ssp585/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmax___ALLRANK01___ssp585.nc']

mf_tasmin = [root_data_dir + 'historical/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmin___ALLRANK01___historical.nc',
             root_data_dir +     'ssp245/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmin___ALLRANK01___ssp245.nc',
             root_data_dir +     'ssp370/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmin___ALLRANK01___ssp370.nc',
             root_data_dir +     'ssp585/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmin___ALLRANK01___ssp585.nc']

mf_pr     = [root_data_dir + 'historical/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYSUM___pr___ALLRANK01___historical.nc',
             root_data_dir +     'ssp245/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYSUM___pr___ALLRANK01___ssp245.nc',
             root_data_dir +     'ssp370/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYSUM___pr___ALLRANK01___ssp370.nc',
             root_data_dir +     'ssp585/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYSUM___pr___ALLRANK01___ssp585.nc']


mf_hist = [root_data_dir + 'historical/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmax___ALLRANK01___historical.nc',
           root_data_dir + 'historical/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmin___ALLRANK01___historical.nc',
           root_data_dir + 'historical/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYSUM___pr___ALLRANK01___historical.nc']


mf_rc45 = [root_data_dir + 'ssp245/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmax___ALLRANK01___ssp245.nc',
           root_data_dir + 'ssp245/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmin___ALLRANK01___ssp245.nc',
           root_data_dir + 'ssp245/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYSUM___pr___ALLRANK01___ssp245.nc']


mf_rc70 = [root_data_dir + 'ssp370/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmax___ALLRANK01___ssp370.nc',
           root_data_dir + 'ssp370/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmin___ALLRANK01___ssp370.nc',
           root_data_dir + 'ssp370/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYSUM___pr___ALLRANK01___ssp370.nc']


mf_rc85 = [root_data_dir + 'ssp585/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmax___ALLRANK01___ssp585.nc',
           root_data_dir + 'ssp585/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN___tasmin___ALLRANK01___ssp585.nc',
           root_data_dir + 'ssp585/LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYSUM___pr___ALLRANK01___ssp585.nc']

#
##########################################################


# In[ ]:





# In[5]:


##########################################################
#
# Pull Historical 
#

ds_tasmax = xr.open_mfdataset(paths      = mf_tasmax,
                              combine    = 'nested',
                              join       = 'outer',
                              concat_dim = scenario).sel(year=my_years) 

ds_tasmin = xr.open_mfdataset(paths      = mf_tasmin,
                              combine    = 'nested',
                              join       = 'outer',
                              concat_dim = scenario).sel(year=my_years) 

ds_pr     = xr.open_mfdataset(paths      = mf_pr,
                              combine    = 'nested',
                              join       = 'outer',
                              concat_dim = scenario).sel(year=my_years) 

tasavg                        = (ds_tasmax["tasmax"]+ds_tasmin["tasmin"])/2
#ds_tasavg.name    = "tasavg"
tasavg.attrs["long_name"]     = "2-m Mean Daily Air Temperature"
tasavg.attrs["description"]   = "2-m Mean Daily Air Temperature"
tasavg.attrs["cell_methods"]  = "time: mean within days  time: mean over months"
tasavg.attrs["standard_name"] = "air_temperature"
tasavg.attrs["units"]         = "degC"
print("Making Tasavg Dataframe")
ds_tasavg = xr.Dataset(data_vars = {"tasavg": tasavg})

print(ds_tasmin)
print("")
print("########################################")
print("")

print(ds_tasmax)
print("")
print("########################################")
print("")

print(ds_tasavg)
print("")
print("########################################")
print("")

print(ds_pr)
print("")
print("########################################")
print("")
#
if Truw:
  ds_tasmax.to_netcdf(path           =  "./tasmax.nc",
               mode           =            'w', 
               format         =      "NETCDF4",
               #engine         =     "h5netcdf",
               encoding       = {"tasmax": {#           "zlib":    True,
                                            #     "complevel" :       7, 
                                                      "dtype": "int16", 
                                               "scale_factor":     0.1,
                                                 "add_offset":     0.0,                                        
                                                 "_FillValue":  -32767}})
  print("tasmax")

  ds_tasmin.to_netcdf(path           =  "./tasmin.nc",
               mode           =            'w', 
               format         =      "NETCDF4",
               engine         =     "h5netcdf",
               encoding       = {"tasmin": {#           "zlib":    True,
                                            #     "complevel" :       7, 
                                                      "dtype": "int16", 
                                               "scale_factor":     0.1,
                                                 "add_offset":     0.0,                                        
                                                 "_FillValue":  -32767}})
  print("tasmin")

  ds_pr.to_netcdf(path           =  "./pr.nc",
               mode           =            'w', 
               format         =      "NETCDF4",
               engine         =     "h5netcdf",
               encoding       = {"pr": {               "zlib":    True,
                                                 "complevel" :       7, 
                                                      "dtype": "int16", 
                                               "scale_factor":     0.1,
                                                 "add_offset":     0.0,                                        
                                                 "_FillValue":  -32767}})
  print("pr")


print("000000")


ds_tasavg.to_netcdf(path           =  "./tasavg.nc",
             mode           =            'w', 
             format         =      "NETCDF4",
             engine         =     "h5netcdf",
             encoding       = {"tasavg": {           "zlib":    True,
                                               "complevel" :       7, 
                                                    "dtype": "int16", 
                                             "scale_factor":     0.1,
                                               "add_offset":     0.0,                                        
                                               "_FillValue":  -32767}})
print("tasavg")




# In[ ]:




