#!/usr/bin/env python
# coding: utf-8

# # LOCA2 30-year Moving Mean Monthly Min Temps Historical Period

# In[ ]:


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

import datetime as datetime

# Loading pandas

import pandas            as pd

import os as os

import subprocess as subprocess

def geo_idx(dd, dd_array):

    geo_idx = (np.abs(dd_array - dd)).argmin()
    return geo_idx

#
##########################################################


# ##  File Control

# In[ ]:





# In[ ]:


##########################################################
#
# File Control
#

Original_File_Prefix = "LOCA2-CONUS-MONTHLY_MEAN"
Final_File_Prefix    = "LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN"

variable    = "tasmin"

tempfile    = "./" + variable + "_tempfile_h.nc"
memberfile  = "./" + variable + "_model_member_h.nc"

local_hdf_string = "export HDF5_USE_FILE_LOCKING=FALSE && "
local_hdf_string = " "

cell_method    = "time: minimum within days  time: mean within months  time: mean over 30 years "

target_rank =  "1"
rank00      = "01"

root_directory        = "/data/DATASETS/LOCA_MACA_Ensembles/LOCA2/LOCA2_CONUS/Climate_CONUS/Monthly"
root_url              = "http://kyrill.ias.sdsmt.edu:8080/thredds/dodsC/LOCA2/Climate_CONUS/Monthly"

loca2_inventory_file  = "/data/DATASETS/LOCA_MACA_Ensembles/LOCA2/LOCA2_CONUS/LOCA2_Model_Member_Available_List.csv"
loca2_complete_file   = "/data/DATASETS/LOCA_MACA_Ensembles/LOCA2/LOCA2_CONUS/LOCA2_Model_Member_Complete_List.csv"

loca2_mask            = "/data/DATASETS/LOCA_MACA_Ensembles/LOCA2/LOCA2_CONUS/LOCA2_MASKS.nc"

#
##########################################################


# In[ ]:


##########################################################
#
# Time Coordinates
#

years_start = np.arange(start =    1950, 
                        stop  =    1986,   dtype=np.float32)
years_end   = np.array(years_start + 29, dtype=np.float32) 

years_middle = np.array(years_start/2 + years_end/2, dtype=np.float32)
years_bounds = np.array([years_start,years_end]).transpose()
n_runnings = len(years_start)

yearall = xr.DataArray(name = "year",
                      data   = years_middle,
                      dims   = {"year": n_runnings},
                      coords = {"year": years_middle},
                      attrs  = {"description": "middle calendar year for 30-yr period",
                                "long_name": "middle calendar year for 30-yr period",
                                "bounds":"year_bnds"})


year_bnds = xr.DataArray(name = "year_bnds",
                      data   =  years_bounds,
                      coords = {"year": years_middle,
                                "bnds": 2},
                      attrs  = {"description": "boundary calendar years for 30-yr period",
                                "long_name": "boundary calendar years for 30-yr period"})


ds_masks = xr.open_dataset(filename_or_obj = loca2_mask)

lon_bnds = ds_masks["lon_bounds"]
lat_bnds = ds_masks["lat_bounds"]
lat_bnds.name = "lat_bnds"
lon_bnds.name = "lon_bnds"

lon = ds_masks["lon"]
lat = ds_masks["lat"]


print(lat)
print(lon)
print(year_bnds)

#
##########################################################


# In[ ]:





# ## Inventories and Lookup Tables
# ### All Potential Ensembles

# In[ ]:


##########################################################
#
# Inventory and Model/Member Lookup Table for all Possible Ensembles
#


loca2_complete_list = pd.read_csv(filepath_or_buffer = loca2_complete_file)

modelsc               = loca2_complete_list[         "Model"].values
membersc              = loca2_complete_list[        "Member"].values
model_membersc        = loca2_complete_list[  "Model_Member"].values

model_member_key0 = np.array(modelsc+"."+membersc, dtype="str")
model_member_key =  np.array(["NULL"], dtype="str")

model_members_to_save = np.append(0,model_membersc)

model_member_key = np.append(model_member_key,model_member_key0)


model_member_key    = xr.DataArray(model_member_key, 
                                   coords={"model_member":model_members_to_save},
                                   name  = "model_member_key",
                                   dims  = ["model_member"],
                                   attrs = {"description" : "Model and Member Label",
                                              "long_name" : "Model and Member Label",
                                               "comment1" :   "LUT Indexing Starts at 0"})

print(model_member_key)

#
##########################################################


# ### All Available Ensembles for Selected Rank

# In[ ]:


##########################################################
#
# Inventory and Model/Member Lookup Table for all Possible Ensembles
#

loca2_ensembles_list = pd.read_csv(filepath_or_buffer = loca2_inventory_file)

loca2_ensembles_list = loca2_ensembles_list.query("Rank == " + target_rank)

models               = loca2_ensembles_list[         "Model"].values
members              = loca2_ensembles_list[        "Member"].values
model_members        = loca2_ensembles_list[  "Model_Member"].values
n_complete_enss      = loca2_ensembles_list["n_complete_ens"].values
historical_invs      = loca2_ensembles_list[    "historical"].values
ssp245_invs          = loca2_ensembles_list[        "ssp245"].values
ssp370_invs          = loca2_ensembles_list[        "ssp370"].values
ssp585_invs          = loca2_ensembles_list[        "ssp585"].values
prec_invs            = loca2_ensembles_list[            "pr"].values
tmax_invs            = loca2_ensembles_list[        "tasmax"].values
tmin_invs            = loca2_ensembles_list[        "tasmin"].values

scenarios            = ["historical", 
                            "ssp245", 
                            "ssp370", 
                            "ssp585"]

print(loca2_ensembles_list)

#
##########################################################


# ## Loop for Results

# In[ ]:


##########################################################
#
# Loop Test
#


for scenario in scenarios[0:1]:
    print("# ################################################")
    print("# ################################################")
    print("# ################################################")
    First = True
    for m in range(len(models)-1):
        print("# ------------------------------------------------")

        model          =          models[m]
        member         =         members[m]
        model_member   =   model_members[m]

        n_complete_ens = n_complete_enss[m]
        prec_inv       =       prec_invs[m]
        tmax_inv       =       tmax_invs[m]
        tmin_inv       =       tmin_invs[m]

        model_member_name = np.array([model + "." + member], dtype="str")
        model_member      = np.array([model_member], dtype="int16").flatten().astype(np.int16)



        model_member = xr.DataArray(data   = model_member.astype(np.int16),
                                    coords = {"model_member":model_member.astype(np.int16)},
                                    name   =  "model_member",
                                    dims   = ["model_member"],
                                    attrs  = {"description" : "Model and Member Code",
                                              "long_name"   : "Model and Member Code",
                                              "code_to_name_lookup_table":  model_member_key.values.tolist(),
                                              "comment1"    : "LUT Indexing Starts at 0"})

        print("# ================================================")

        # print(m,scenario,loca2_ensembles_list.iloc[m])

        inventory = loca2_ensembles_list.iloc[m].loc[scenario]
        print("# " + str(model_member[0].values).zfill(3) + " " + model + " " + member + " " + scenario + " " + inventory)

        if (inventory != "---"):


            if ("N" in inventory):



                print("#  . . . . . . . . . . . . . . . . . . . . . . . .")

                hist_file         = root_directory                            +  "/"  + \
                                    "historical"                              +  "/"  + \
                                    Original_File_Prefix                      + "___" + \
                                    variable                                  + "___" + \
                                    models[m]                                 +  "."  + \
                                    members[m]                                + "___" + \
                                    "historical"                              + ".nc"  



                combined_file     = root_directory                            +  "/"  + \
                                    Final_File_Prefix                         + "___" + \
                                    variable                                  + "___" + \
                                    str(model_member[0].values).zfill(3)      + "___" + \
                                    models[m]                                 +  "."  + \
                                    members[m]                                + "___" + \
                                    scenario                                  + ".nc"  

                combined_wc_files = root_directory                            +  "/"  + \
                                    Final_File_Prefix                         + "___" + \
                                    variable                                  + "___" + \
                                    "???"                                     + "___" + \
                                    "*"                                       + "___" + \
                                    scenario                                  + ".nc" 

                final_merged_file = root_directory                            +  "/"  + \
                                    scenario                                  +  "/"  + \
                                    Final_File_Prefix                         + "___" + \
                                    variable                                  + "___" + \
                                    "ALLRANK"+ rank00                         + "___" + \
                                    scenario                                  + ".nc" 


                print("    - hist_file: " + hist_file)
                print("    - comb_file: " + combined_file)
                print("    - comw_file: " + combined_wc_files)
                print("    - finl_file: " + final_merged_file)


                ds            = xr.open_dataset(filename_or_obj = hist_file)
                time_hist_max = ds["time"].values.max()
                time_hist_n   = ds[variable].values.shape

                print("#    Max_Orig_Time = " + str(time_hist_max) + "      " + str(time_hist_n) )

                if (First):
                    model_member_array = np.array(model_members[m], dtype = "int16")
                    First              = False
                else:
                    model_member_array = np.append(model_member_array, model_members[m]).flatten().astype(np.int16)



                cdo_cat_command = "cdo --no_history -f nc4 -z zip_9  mergetime "

                command_aggregate = "cp -frv " + hist_file + " " + tempfile
                subprocess.run(["rm -fr " + tempfile],                             shell = True, check = True)
                subprocess.run([command_aggregate],                                shell = True, check = True)
                subprocess.run(["ncatted -Oh -a bounds,lon,d,,      " + tempfile], shell = True, check = True)
                subprocess.run(["ncatted -Oh -a bounds,lat,d,,      " + tempfile], shell = True, check = True)
                subprocess.run(["ncatted -Oh -a _FillValue,lon,d,,  " + tempfile], shell = True, check = True)
                subprocess.run(["ncatted -Oh -a _FillValue,lat,d,,  " + tempfile], shell = True, check = True)
                subprocess.run(["ncatted -Oh -a _FillValue,time,d,, " + tempfile], shell = True, check = True)

                ds              = xr.open_dataset(filename_or_obj = tempfile)
                tasmin0         = ds[variable]

                time_merged_max = ds["time"].values.max()
                time_merged_n   = ds[variable].values.shape

                print("#   Max_Merge_Time = " + str(time_merged_max) + "     " + str(time_merged_n) )      

                subprocess.run(["rm -fr " + tempfile], shell = True, check = True)



                time = ds["time"]
                lon = ds["lon"]
                lat = ds["lat"]
                nt = time.shape[0]
                nm = 12
                ny = nt/12


                start_year = time.dt.year.min().values
                end_year   = time.dt.year.max().values
                year  = np.arange(start = start_year,
                                  stop  = end_year+1,
                                  dtype = np.int16)
                month = np.arange(start = 1,
                                  stop  = 12+1,
                                  dtype = np.int16)

                yeardv = xr.DataArray(name = "year",
                                      data   =  year,
                                    dims   = "year",
                                    coords = {"year": year},
                                    attrs  = {"description": "beginning period calendar year",
                                                "long_name": "beginning period calendar year"})



                monthdv = xr.DataArray(data   =  month,
                                      dims    = "month",
                                      coords  = {"month": month},
                                      attrs   = {"description": "calendar month",
                                                   "long_name": "calendar month"})

                print("Start Reindexing",os.system("date"))

                multiindex_ds = ds.assign_coords(month = monthdv,
                                                 year  = yeardv).  \
                                   stack(time2d=("year",
                                                 "month")).        \
                                   reset_index("time", drop=True). \
                                   rename(time="time2d").          \
                                   unstack("time2d").drop_vars("time_bnds")

                print("Start Rolling Mean",os.system("date"))


                for i in range(n_runnings):
                    if ((i % 10) == 0):
                        print("   --- Processing ",years_start[i] , "to", years_end[i])

                    year_co = xr.DataArray(name = "year",
                                           data   =np.array([years_start[i]/2+years_end[i]/2], dtype=np.float32),
                                           dims   = "year",
                                           coords = {"year": np.array([years_start[i]], dtype=np.int16)},
                                           attrs  = {"description": "beginning period calendar year",
                                                     "long_name": "beginning period calendar year"})

                    if (i == 0):
                        running_var = multiindex_ds[variable].sel(year = slice(years_start[i],years_end[i])).mean(dim="year", keep_attrs = True).transpose("month","lat","lon").expand_dims(dim={"model_member" : 1,"year":1}) 
                        running_var.coords["model_member"]=model_member
                        running_var.coords["year"] = year_co
                    else:
                        temp_30 = multiindex_ds[variable].sel(year = slice(years_start[i],years_end[i])).mean(dim="year").transpose("month","lat","lon").expand_dims(dim={"model_member" : 1,"year":1}) 
                        temp_30.coords["model_member"]=model_member
                        temp_30.coords["year"] = year_co
                        running_var = xr.concat([running_var, temp_30], dim = "year")
                        del temp_30

                print("Finished Rolling Mean",os.system("date"))


                outdata = xr.Dataset(data_vars = {"model_member" : model_member.astype(np.int16),
                                                  "year"         : yearall,
                                                #  "year_bnds"    : year_bnds,
                                                  "month"        : monthdv,
                                                  "lat"          : lat,
                                                  #"lat_bnds"     : lat_bnds,
                                                  "lon"          : lon,
                                                  #"lon_bnds"     : lon_bnds,
                                                   variable      : running_var},
                                     attrs     = {"scenario"     : scenario})

                running_var.attrs["cell_methods"] = cell_method



                outdata.to_netcdf(path           =  combined_file, 
                                  mode           =            'w', 
                                  format         =      "NETCDF4",
                                  engine         =     "h5netcdf", #
                                  unlimited_dims = "model_member",
                                  encoding       = {variable: {         "zlib":    True,
                                                                  "complevel" :       7, 
                                                                       "dtype": "int16", 
                                                                "scale_factor":     0.1,
                                                                  "add_offset":     0.0,                                        
                                                                  "_FillValue":  -32767}})

                print("Writing NetCDF Climate File",os.system("date"))




            # end check on available variable
        # end check on on available member

    #end loop on model

    print("# = = = = = = = = = = = = = = = = = = = = = = = = ") 

    model_member = xr.DataArray(data   = model_member_array.astype(np.int16),
                                name   =  "model_member",
                                dims   = ["model_member"],
                                attrs  = {"description" : "Model and Member Code",
                                          "long_name"   : "Model and Member Code",
                                          "code_to_name_lookup_table":  model_member_key.values.tolist(),
                                          "comment1"    : "LUT Indexing Starts at 0"})

    model_member_ds = xr.Dataset(data_vars = {"model_member" : model_member})

    model_member_ds.to_netcdf(path            =   memberfile, 
                               mode           =           'w', 
                               format         =     "NETCDF4",
                               engine         =    "h5netcdf", #
                               unlimited_dims = "model_member")  

    cdo_cat_command = " cdo --no_history -f nc4 -z zip_9 cat "
    nco_cat_command = " ncrcat --4 --hst --dfl_lvl 9  "
    command_aggregate = nco_cat_command +combined_wc_files + " " + final_merged_file    
    print("# Final Aggregation for "+scenario)
    subprocess.run(["rm -fr " + final_merged_file + " " + tempfile + " 2_" + tempfile], 
                   shell = True, 
                   check = True)    
    subprocess.run([local_hdf_string + command_aggregate], 
                   shell = True, 
                   check = True)
    print("# Files Concatenated")
    subprocess.run(["rm -frv " + combined_wc_files], 
                   shell = True, 
                   check = True)
    print("# Files Cleaned")


# end loop on scenario

print("# ================================================")
print("end processing")


# In[ ]:





# In[ ]:





# In[ ]:




