#!/usr/bin/env python
# coding: utf-8

# # LOCA2 30-year Moving Mean Annual Max Temps.

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


##########################################################
#
# File Control
#

Original_File_Prefix = "LOCA2-CONUS-MONTHLY_MEAN"
Final_File_Prefix    = "LOCA2-CONUS-ANNUAL30YRUNMEAN_MONTHLYMEAN"

variable    = "tasmax"

tempfile    = "./" + variable + "_tempfile.nc"
memberfile  = "./" + variable + "_model_member.nc"

local_hdf_string = "export HDF5_USE_FILE_LOCKING=FALSE && "
local_hdf_string = " "

cell_method    = "time: maximum within days  time: mean within years  time: mean over 30 years "
cell_methodsdv = "time: mean over months   time: stdev over 30 years "

target_rank =  "1"
rank00      = "01"

root_directory        = "/data/DATASETS/LOCA_MACA_Ensembles/LOCA2/LOCA2_CONUS/Climate_CONUS/Monthly"
root_url              = "http://kyrill.ias.sdsmt.edu:8080/thredds/dodsC/LOCA2/Climate_CONUS/Monthly"

loca2_inventory_file  = "/data/DATASETS/LOCA_MACA_Ensembles/LOCA2/LOCA2_CONUS/Original_CONUS/LOCA2_Model_Member_Available_List.csv"
loca2_complete_file   = "/data/DATASETS/LOCA_MACA_Ensembles/LOCA2/LOCA2_CONUS/Original_CONUS/LOCA2_Model_Member_Complete_List.csv"

#
##########################################################


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

    
for scenario in scenarios[1:]:
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

                futr_file         = root_directory                            +  "/"  + \
                                    scenario                                  +  "/"  + \
                                    Original_File_Prefix                      + "___" + \
                                    variable                                  + "___" + \
                                    models[m]                                 +  "."  + \
                                    members[m]                                + "___" + \
                                    scenario                                  + ".nc"  

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
        
                ds            = xr.open_dataset(filename_or_obj = futr_file)
                time_futr_max = ds["time"].values.max()
                time_futr_n   = ds[variable].values.shape
                
                print("#    Max_Orig_Time = " + str(time_futr_max) + "      " + str(time_futr_n) )
                
                if (First):
                    model_member_array = np.array(model_members[m], dtype = "int16")
                    First              = False
                else:
                    model_member_array = np.append(model_member_array, model_members[m]).flatten().astype(np.int16)
                    
                    

                cdo_cat_command = "cdo --no_history -f nc4 -z zip_9  mergetime "

                command_aggregate = cdo_cat_command + hist_file + " " + futr_file + " " + tempfile
                subprocess.run(["rm -fr " + tempfile],                             shell = True, check = True)
                subprocess.run([command_aggregate],                                shell = True, check = True)
                subprocess.run(["ncatted -Oh -a bounds,time,d,,     " + tempfile], shell = True, check = True)
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
                
                yeardv = xr.DataArray(data   =  year,
                                    dims   = "year",
                                    coords = {"year": year},
                                    attrs  = {"description": "calendar year",
                                                "long_name": "calendar year"})
                
                monthdv = xr.DataArray(data   =  month,
                                      dims    = "month",
                                      coords  = {"month": month},
                                      attrs   = {"description": "calendar month",
                                                   "long_name": "calendar month"})
                
                print("Start Reindexing",os.system("date"))
                
                multiindex_ds = ds.assign_coords(month =monthdv,
                                                 year  =yeardv
                ).stack(
                    time2d=("year","month")
                ).reset_index(
                    "time", drop=True
                ).rename(
                    time="time2d"
                ).unstack("time2d")

                print("Start Rolling Mean",os.system("date"))

                
                rolling_monthly = multiindex_ds[variable].rolling(year   =   30,
                                                                  center = True).mean().dropna(dim = "year",
                                                                                               how =  "all") 
                rolling_monthly.expand_dims(dim={"model_member" : 1}) 
                rolling_monthly.attrs["cell_methods"] = cell_method
                #del rolling_monthly.attrs["coordinates"]
                print("Finished Rolling Mean",os.system("date"))


                outdata = xr.Dataset(data_vars = {variable       :               rolling_monthly,
                                                  "lat"          :                           lat,
                                                  "lon"          :                           lon,
                                                  "model_member" : model_member.astype(np.int16)},
                                     attrs     = {"scenario"     :                      scenario})

                


                outdata.to_netcdf(path           =  combined_file, 
                                  mode           =            'w', 
                                  format         =      "NETCDF4",
                                  engine         =     "h5netcdf", #
                                  unlimited_dims = "model_member",
                                  encoding       = {variable: {        "dtype": "int16", 
                                                                "scale_factor":     0.1,
                                                                  "add_offset":     0.0,                                        
                                                                  "_FillValue":  -32767}})

                print("Writing NetCDF Climate File",os.system("date"))
                
                subprocess.run([local_hdf_string+" ncatted -Oh -a _FillValue,lon,d,, " + combined_file], 
                               shell = True, 
                               check = True)
                subprocess.run([local_hdf_string+" ncatted -Oh -a _FillValue,lat,d,, " + combined_file], 
                               shell = True, 
                               check = True)

                subprocess.run([local_hdf_string + " ncpdq -h -a year,month,lat,lon " + combined_file + " " + combined_file+".swapped.nc"],
                                shell = True, 
                                check = True)      
                print("Correcting all Dimensions",os.system("date"))

                subprocess.run(["mv -v " + combined_file + ".swapped.nc " + combined_file],
                                shell = True,  
                                check = True)                               
                ds               = xr.open_dataset(filename_or_obj = combined_file)
                time_running_max = ds["year"].values.max()
                time_running_n   = ds[variable].values.shape
                print("# Max_Running_Time = " + str(time_running_max) + "  " + str(time_running_n) )                                 
                
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
    nco_cat_command = " ncecat -h -M -u model_member "
    command_aggregate = cdo_cat_command +combined_wc_files + " " + tempfile    
    print("# Final Aggregation")
    subprocess.run(["rm -fr " + final_merged_file + " " + tempfile + " 2_" + tempfile], 
                   shell = True, 
                   check = True)    
    subprocess.run([local_hdf_string + command_aggregate], 
                   shell = True, 
                   check = True)
    print("# Files Concatenated")

    subprocess.run([local_hdf_string+" ncks -h -C -O -x -v model_member " + tempfile + " " + final_merged_file], 
                   shell = True, 
                   check = True)
    print("# dimension dropped")
    subprocess.run([local_hdf_string+" ncks -h -A "+ memberfile + " " + final_merged_file], 
                   shell = True, 
                   check = True)
    ds               = xr.open_dataset(filename_or_obj = final_merged_file)
    time_running_max = ds["year"].values.max()
    time_running_n   = ds[variable].values.shape
    print("#     Max_Ens_Time = " + str(time_running_max) + " " + str(time_running_n) )  
    print("# dimension swapped")
    subprocess.run(["rm -fr  " + tempfile + " 2_" + tempfile + " " + combined_wc_files], 
                   shell = True, 
                   check = True) 

# end loop on scenario
                
print("# ================================================")
print("end processing")


# In[ ]:





# 

# In[ ]:




