# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 12:06:56 2021

@author: judithg
"""
import os
import netCDF4
from netCDF4 import num2date
import numpy as np
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
    
    for i in range(dt,length_dt,dt):
        print(i)
        nc_fnameX = nc_fname + '_' +  "{0:0>3}".format(i) +'.nc'
        nc = netCDF4.Dataset(nc_fnameX)
        
        tt = np.append(tt,num2date(np.array(nc['time']), units=nc['time'].units, calendar = nc['time'].calendar),axis=0)
        Uwind = np.append(Uwind,np.array(nc['u-component_of_wind_height_above_ground']),axis=0)
        Vwind = np.append(Vwind,np.array(nc['v-component_of_wind_height_above_ground']),axis=0)
        print(nc_fnameX)
    Uwind = np.squeeze(Uwind)
    Vwind = np.squeeze(Vwind)
    return(ln,lt,tt,Uwind,Vwind)
    
def write_wind_to_ascii(U,V,t,lnn,ltt,out_file):
    file = open(out_file,"w")
    for t in range(len(t)):
        for m in range(len(ltt)):
            for n in range(len(lnn)):
                file.write(str("%.2f" % U[t,m,n]) + "\t")
            file.write("\n")
        for m in range(len(ltt)):
            for n in range(len(lnn)):
                file.write(str("%.2f" % V[t,m,n]) + "\t")
            file.write("\n")               
    file.close()
    print(out_file + ' created')
    return()

def make_wind_file(fl_name,U,V,t,lnn,ltt):
    ## interpolate wind data onto time series with equal spacing
    #t_old = t
    #U_new,V_new,t_new = interp_wnd_2_lin_timeseries(U,V,t_old,0.25)
    ## save wind data to ascii file
    out_file = fl_name + '.dat'
    #print('creating ' + out_file)
    write_wind_to_ascii(U,V,t,lnn,ltt,out_file)
    return()

###############################################################################
def make_winds(now):
    
    
    wnd_nc = '../tmp/wind_tmp'
    hwnd_nc = '../tmp/h_wind_tmp'
    hhwnd_nc = '../tmp/hh_wind_tmp'
    wind_dt = 3
    wind_length = 181
    hwind_length = 22
    
    # now = datetime.datetime.utcnow()- timedelta(0.125)# NCEP needs at least 3 hours to upload their forecast from UTC00
    out_name='../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") + now.strftime("%H")
    print('Gnerating wind forcing in: ' + out_name)
    try:
        os.mkdir(out_name)
    except OSError as error:
        print(error)    
    out_name='../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H")  + '/winds'
    
    ln,lt,tt,Uwind,Vwind = read_netcdf(wnd_nc,wind_dt,wind_length)
    ln,lt,htt,hUwind,hVwind = read_netcdf(hwnd_nc,wind_dt,hwind_length)
    ln,lt,hhtt,hhUwind,hhVwind = read_netcdf(hhwnd_nc,wind_dt,hwind_length)
    
    # TT=np.concatenate((htt,tt), axis=0) 
    # time_str=[TT[i].strftime('%Y-%m-%d %H:%M') for i in range(0,len(tt))]
    
    ttt_size = (len(tt)+len(htt)+len(hhtt))*3
    ttt = np.arange(0, ttt_size, 3.0)
    
    tUwind = np.concatenate((hhUwind,hUwind,Uwind), axis=0)
    tVwind = np.concatenate((hhVwind,hVwind,Vwind), axis=0)
    
    
    make_wind_file(out_name,tUwind,tVwind,ttt,ln,lt)
    




