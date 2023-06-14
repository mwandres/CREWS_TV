import scipy.io as sio
import numpy as np
from matplotlib import pyplot as plt
import datetime
import pandas as pd
import netCDF4 as nc
from datetime import timedelta

def readvars_nc(nc_fname,varname):
    ncc = nc.Dataset(nc_fname,'r') 

    if varname == 'time':
        time_var = ncc.variables[varname] 
        time_or = nc.num2date(ncc[varname],time_var.units, time_var.calendar)
        time_str=[time_or[i].strftime('%Y-%m-%d %H:%M') for i in range(0,len(time_or))]
        var = [datetime.datetime.strptime(time_str[i],'%Y-%m-%d %H:%M') for i in range(0,len(time_str))]
    else:
        var = np.array(ncc[varname])
        
    ncc.close()
    return(var)

def generate_tides(name,start_date,folder_tmp):
    nc_fname = folder_tmp  + 'tide_hourly_'+name+'.nc'
    time_tide = readvars_nc(nc_fname,'time')
    tide = readvars_nc(nc_fname,'tide')/100

    nc_fname = folder_tmp  + 'tide_minute_'+name+'.nc'
    time_tide_min = readvars_nc(nc_fname,'time')
    tide_min = readvars_nc(nc_fname,'tide')/100
    time_w = readvars_nc(nc_fname,'time')

    xor = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_w]
    xnew = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_tide_min]
    diff = np.zeros_like(time_tide_min, dtype=np.float64)
    #-- forward differentiation for starting point
    diff[0] = tide_min[1] - tide_min[0]
    #-- backward differentiation for end point
    diff[-1] = tide_min[-1] - tide_min[-2]
    #-- centered differentiation for all others
    diff[1:-1] = (tide_min[2:] - tide_min[0:-2])/2.0
    htindex, = np.nonzero((np.sign(diff[0:-1]) >= 0) & (np.sign(diff[1:]) < 0))
    ltindex, = np.nonzero((np.sign(diff[0:-1]) <= 0) & (np.sign(diff[1:]) > 0))
    next_hour = time_tide_min[0].replace(microsecond=0, second=0, minute=0) + timedelta(hours=1)
    minutes_diff = (next_hour - time_tide_min[0]).total_seconds() / 60.0

    tidal_time = time_tide_min[int(minutes_diff):]
    tide_height = tide_min[int(minutes_diff):]
    df = pd.DataFrame({'datetime':tidal_time, 'height':tide_height})
    df['datetime'] = pd.to_datetime(df['datetime'])
    mask = df['datetime'] >= start_date
    df2 = df.loc[mask]
    df2.to_csv('tide_data.csv', index=None, header=True)
    return(time_tide_min, tide_min, htindex, ltindex)