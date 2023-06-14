#test
import netCDF4
from netCDF4 import num2date
import datetime as dt 
import calendar
import numpy as np
import xarray as xr
from scipy.interpolate import interp1d
import pandas as pd

def generate_gfs(target_lon, target_lat,nc_fname):
    dt = 0
    length_dt = 181
    visarr=[]
    temparr=[]
    gustarr=[]
    timea=[]
    for i in range(dt,length_dt,3):
        nc_fnameX = nc_fname + '_' +  "{0:0>3}".format(i) +'.nc'
        nc = netCDF4.Dataset(nc_fnameX)
        ds = xr.open_dataset(nc_fnameX)
        dsloc = ds.sel(lon=target_lon,lat=target_lat,method='nearest')
        VIS =  np.array(dsloc['Visibility_surface'][0])
        Temperature_surface =  np.array(dsloc['Temperature_surface'][0])
        celcius = Temperature_surface - 273.15
        Wind_speed_gust_surface =  np.array(dsloc['Wind_speed_gust_surface'][0])
        timer = num2date(np.array(nc['time']), units=nc['time'].units, calendar = nc['time'].calendar)
        timea.append(str(timer[0]))
        visarr.append(VIS)
        temparr.append(celcius)
        gustarr.append(Wind_speed_gust_surface)


    df = pd.DataFrame({'datetime':timea, 'vis':visarr, 'sst':temparr, 'wind_gust':gustarr})
    df.to_csv('gfs_3hourly.csv', index=None, header=True)
    return None

def get_start_date():
    df = pd.read_csv('gfs_3hourly.csv')
    df['datetime'] = pd.to_datetime(df['datetime'])
    return df.iloc[0]['datetime']

def toTimestamp(d): 
    return calendar.timegm(d.timetuple()) 

def interpolate_gfs():
    #read files 
    data_hourly = pd.read_csv('dataset.csv') 
    data_3hourly = pd.read_csv('gfs_3hourly.csv')

    dates_hourly_list = [dt.datetime.strptime(date, "%Y-%m-%d %H:%M:%S") for date in data_hourly.DateTime] 
    dates_3hourly_list = [dt.datetime.strptime(date, "%Y-%m-%d %H:%M:%S") for date in data_3hourly.datetime] 
    sst_3hrly = data_3hourly.sst.values
    vis_3hrly = data_3hourly.vis.values
    gust_3hrly = data_3hourly.wind_gust.values

    # create time arrays 
    t1 = np.zeros(len(dates_hourly_list)) 
    for i in range(len(t1)):
        t1[i] = toTimestamp(dates_hourly_list[i])

    t2 = np.zeros(len(dates_3hourly_list)) 
    for i in range(len(t2)):
        t2[i] = toTimestamp(dates_3hourly_list[i])
        
    # interpolate sst
    f = interp1d(t2, sst_3hrly,bounds_error=False, kind='nearest') 
    sst_1hrly = f(t1)

    f = interp1d(t2, vis_3hrly,bounds_error=False, kind='nearest') 
    vis_1hrly = f(t1)

    f = interp1d(t2, gust_3hrly,bounds_error=False) 
    gust_1hrly = f(t1)

    dfx = pd.DataFrame({'datetime':dates_hourly_list, 'vis':np.round(vis_1hrly,decimals = 2), 'sst':np.round(sst_1hrly,decimals = 2), 'wind_gust':np.round(gust_1hrly,decimals = 2)})
    dfx.to_csv('gfs_3hourly_interpolated.csv', index=None, header=True)