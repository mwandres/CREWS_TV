# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 03:29:44 2022

@author: moritzw
"""
import xarray as xr
import datetime as dt
import wavespectra
import numpy as np

def load_wave_spectra(now,name_spec):
    fl_name_spec_1 = '../runs/'+now.strftime('%Y%m%d%H')+'/'+ name_spec +'.spec'
    sp = wavespectra.read_swan(fl_name_spec_1)
    return(sp)

def create_superpoint(sp1,sp2,sp3,sp4):
    swan_time = sp1.efth.time.data
    swan_dir = sp1.efth.dir.data
    swan_freq = sp1.efth.freq.data
    efth1 = np.squeeze(sp1.efth.data)
    efth2 = np.squeeze(sp2.efth.data)
    efth3 = np.squeeze(sp3.efth.data)
    efth4 = np.squeeze(sp4.efth.data)
    efth_super = np.zeros(np.shape(efth1))
    
    for t in range(len(swan_time)):
        for f in range(len(swan_freq)):
            for d in range(len(swan_dir)):
                tmp = (efth1[t,f,d],efth2[t,f,d],efth3[t,f,d],efth4[t,f,d])
                efth_super[t,f,d] = np.max(tmp)
                
    return(swan_time,swan_dir,swan_freq,efth_super)

def load_swan_output(now,fl_name):
    swan_table = np.loadtxt(
        '../runs/'+now.strftime('%Y%m%d%H')+'/'+ fl_name,
        skiprows = 7)
    return(swan_table)

def extract_avg_wind_from_tables(
        swan_table_1,swan_table_2,swan_table_3,swan_table_4):
    
    wind_x_all = np.array(
        [swan_table_1[:,7],swan_table_2[:,7],
         swan_table_3[:,7],swan_table_4[:,7]]
        ).T
    
    wind_y_all = np.array(
        [swan_table_1[:,8],swan_table_2[:,8],
         swan_table_3[:,8],swan_table_4[:,8]]
        ).T
    
    wind_x = np.mean(wind_x_all,axis=1)
    wind_y = np.mean(wind_y_all,axis=1)
    return(wind_x,wind_y)

def save_data_as_nc(nc_out_name,
        efth_super,wind_x,wind_y,swan_dir,swan_freq,swan_time):
    df1 = xr.DataArray(efth_super.T, coords=[('dir', swan_dir), ('freq', swan_freq),('time',swan_time)])
    df1.name = 'efth'
    df2 = xr.DataArray(wind_x,coords=[('time',swan_time)])
    df2.name = 'wind_x'
    df3 = xr.DataArray(wind_y,coords=[('time',swan_time)])
    df3.name = 'wind_y'
    df = xr.merge([df1,df2,df3])
    df.to_netcdf(nc_out_name, engine='netcdf4')
    return()

def create_superpoint_nc(now):
    sp1 = load_wave_spectra(now,'Pto_1')
    sp2 = load_wave_spectra(now,'Pto_2')
    sp3 = load_wave_spectra(now,'Pto_3')
    sp4 = load_wave_spectra(now,'Pto_4')
    
    swan_time,swan_dir,swan_freq,efth_super = create_superpoint(sp1,sp2,sp3,sp4)
    
    swan_table_1 = load_swan_output(now,'Pto_1.tab')
    swan_table_2 = load_swan_output(now,'Pto_2.tab')
    swan_table_3 = load_swan_output(now,'Pto_3.tab')
    swan_table_4 = load_swan_output(now,'Pto_4.tab')
    
    wind_x,wind_y = extract_avg_wind_from_tables(
            swan_table_1,swan_table_2,swan_table_3,swan_table_4)
        
    nc_out_name = '../Funafuti_lagoon/Data/Superpoint_frcst.nc'
    save_data_as_nc(nc_out_name,efth_super,wind_x,wind_y,swan_dir,swan_freq,swan_time)
    return()

# now = dt.datetime(2022,1,6,0)
#create_superpoint_nc(now)    
