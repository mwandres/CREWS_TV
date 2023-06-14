# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 14:43:53 2023

@author: moritzw
"""
import os
import netCDF4
from netCDF4 import num2date
import numpy as np
from matplotlib import pyplot as plt
# import datetime
# from datetime import timedelta

def read_netcdf(nc_fname,dt,length_dt):
    nc_fnameX = nc_fname + '_' + "{0:0>3}".format(0) +'.nc'
    nc = netCDF4.Dataset(nc_fnameX)
    
    ln = np.array(nc['lon'])
    lt = np.array(nc['lat'])
    tt = num2date(np.array(nc['time']), units=nc['time'].units, calendar = nc['time'].calendar)
    Uwind = np.array(nc['u-component_of_wind_height_above_ground'])
    Vwind = np.array(nc['v-component_of_wind_height_above_ground'])
    MSLP = np.array(nc['Pressure_reduced_to_MSL_msl'])
    
    for i in range(dt,length_dt,dt):
        print(i)
        nc_fnameX = nc_fname + '_' +  "{0:0>3}".format(i) +'.nc'
        nc = netCDF4.Dataset(nc_fnameX)
        
        tt = np.append(tt,num2date(np.array(nc['time']), units=nc['time'].units, calendar = nc['time'].calendar),axis=0)
        Uwind = np.append(Uwind,np.array(nc['u-component_of_wind_height_above_ground']),axis=0)
        Vwind = np.append(Vwind,np.array(nc['v-component_of_wind_height_above_ground']),axis=0)
        MSLP =  np.append(MSLP,np.array(nc['Pressure_reduced_to_MSL_msl']),axis=0)

        print(nc_fnameX)
    Uwind = np.squeeze(Uwind)
    Vwind = np.squeeze(Vwind)
    MSLP = np.squeeze(MSLP)
    return(ln,lt,tt,Uwind,Vwind,MSLP)

def plot_winds_and_MSLP(now):
    wnd_nc = '../Regional_tmp/out'
    wind_dt = 12
    wind_length = 181

    ln,lt,tt,Uwind,Vwind,MSLP = read_netcdf(wnd_nc,wind_dt,wind_length)
    
    WndSpd = np.sqrt(Uwind**2+Vwind**2)
    
    xllg = min(ln)
    yllg = min(lt)
    xurg = max(ln)
    yurg = max(lt)
    grid_x = np.linspace(xllg,xurg,20)
    grid_y = np.linspace(yllg,yurg,20)
    lgrid_x, lgrid_y = np.mgrid[xllg:xurg:20j,yllg:yurg:20j]

    ln_x, lt_y = np.meshgrid(ln,lt)
    from scipy.interpolate import interp2d
    for i in range(len(tt)):
        fU = interp2d(ln, lt, Uwind[i,:,:], kind='cubic')
        UU = fU(grid_x, grid_y)
        fV = interp2d(ln, lt, Vwind[i,:,:], kind='cubic')
        VV = fV(grid_x, grid_y)
        
        levels = np.arange(0.,40, 0.1)
        levels_MSLP = np.arange(970.,1030, 10)
        
        
        fig, ax = plt.subplots(figsize=(8,8))
        tcf = ax.contourf(ln,lt,WndSpd[i,:,:],levels=levels,cmap='gist_ncar')
        cbar = fig.colorbar(tcf, ax=ax)
        cbar.set_label('Wind Speed (m/s)',fontsize=16)
        
        ax.quiver(lgrid_x,lgrid_y,UU,VV)
        CS = ax.contour(ln,lt,MSLP[i,:,:]/100,levels=levels_MSLP,colors='w')
        ax.clabel(CS, inline=True, fontsize=14)
        ax.set_xlabel('Longitude (degrees)', fontsize=16)
        ax.set_ylabel('Latitude (degrees)', fontsize=16)
        ax.set_title(str(tt[i]), fontsize=16)
        
        fig_out = '../Regional_Output/out'
        fileprint=fig_out + '_' + "{0:0>3}".format(i) 
        plt.savefig(fileprint)
        plt.close(fig)
    

