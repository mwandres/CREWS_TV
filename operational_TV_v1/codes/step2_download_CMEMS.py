# -*- coding: utf-8 -*-
"""
Created on Thu May  6 17:47:05 2021

@author: judithg
"""
import os
import datetime as dt
import netCDF4 as nc
import numpy as np
import numpy.matlib
import pandas as pd

def read_sla_csv(fpath):
    df = pd.read_csv(fpath)
    atoll = df['atoll'].tolist()
    sla_lon = df['lon'].tolist()
    sla_lat = df['lat'].tolist()
    sla = df['sla_offset'].tolist()
    return(atoll,sla_lon,sla_lat,sla)


###############################################################################
def download_CNEMS(now):  
    print('Downloading Sea Level Anomaly data from http://nrt.cmems-du.eu')
    folder_tmp ='../tmp/' 
    now = now - dt.timedelta(2)# NCEP needs at least 3 hours to upload their forecast from UTC00
    # now = now.replace(hour=0,minute=0,second=0,microsecond=0)
    then = now+dt.timedelta(days=9.5,hours=1)
    
    sla_fpath = '../extras/Tuvalu_island_gauge_locations.csv'
    atoll_list, sla_lon_list, sla_lat_list,sla_offset_list = read_sla_csv(sla_fpath)

    for atoll,sla_lon,sla_lat,sla_off in zip(atoll_list,sla_lon_list,sla_lat_list,sla_offset_list):
        print("\n",atoll,sla_lon,sla_lat)
        
        os.system('python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu\
                              --service-id GLOBAL_ANALYSISFORECAST_PHY_001_024-TDS --product-id cmems_mod_glo_phy_anfc_0.083deg_PT1H-m\
                              --longitude-min %(lon_min)s --longitude-max %(lon_max)s --latitude-min %(lat_min)s --latitude-max %(lat_max)s\
                              --date-min '"'%(ini)s'"' --date-max '"'%(end)s'"'\
                              --depth-min 0.493 --depth-max 0.4942 --variable thetao --variable uo --variable vo --variable zos\
                              --out-dir %(folder)s --out-name sla_%(island)s.nc\
                              --user mwandres --pwd nZp9d@zwPVi73hn' %{"lon_min":(sla_lon-0.001),"lon_max":(sla_lon+0.001),"lat_min":(sla_lat-0.001),"lat_max":(sla_lat+0.001),"ini":now.strftime('%Y-%m-%d %H:%M:%S') ,"end":then.strftime('%Y-%m-%d %H:%M:%S'),"folder":folder_tmp,"island":atoll})
                              
         
                              
        nc_fname = folder_tmp  + ('sla_%s.nc' % atoll)
        ncc = nc.Dataset(nc_fname)    
        time_var = ncc.variables['time']  
        time_or=nc.num2date(time_var[:],time_var.units, time_var.calendar)
        time_str=[time_or[i].strftime('%Y-%b-%d %H:%M') for i in range(0,len(time_or))]
        Time=[dt.datetime.strptime(time_str[i],'%Y-%b-%d %H:%M') for i in range(0,len(time_str))]
        sla = (np.array(ncc['zos'])- sla_off)*100# This offset comes from comparing DUACS altimetry with the model, it has been aplied to transform Sea Surface Height to Sea Level Anomaly (height from the MSL) [1993-2019 monthly mean]
        sla=sla.squeeze()
                
        fn = folder_tmp  + ('sla_hourly_%s.nc' % atoll)
        try:
            os.remove(fn)
        except:
            print('The system cannot find the file specified, creating netcdf file')
        ds = nc.Dataset(fn, 'w', format='NETCDF4')
        time = ds.createDimension('time',None)
        times = ds.createVariable('time', 'f8',('time',))
        times.units='hours since 1950-01-01 00:00:00'
        times.calendar='gregorian'
        SLA = ds.createVariable('SLA','f4', dimensions=('time',))
        SLA.units = 'cm'
        SLA[:]=sla
        times[:]=[nc.date2num(x,units=times.units,calendar=times.calendar)-0.5 for x in Time]
        ds.close()
        
        
        print('sea level anomaly stored as ' + fn)