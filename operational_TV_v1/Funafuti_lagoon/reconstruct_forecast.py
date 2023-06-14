# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 23:39:31 2022

@author: moritzw
"""
#from create_superpoint import create_superpoint_nc
#import datetime as dt
import os
import xarray as xr
import numpy as np
import pickle as pk
from scipy import interpolate as sio
import datetime as dt
from create_superpoint import create_superpoint_nc
import netCDF4 as nc
import matplotlib.dates as mdates
from matplotlib import pyplot as plt
from math import cos, asin, sqrt
import pandas as pd
import numpy.matlib


def determine_clusters(now,path_data):
    num_clusters=2000

    sp_regular=xr.open_dataset(os.path.join(path_data,'Superpoint_frcst.nc'))
    sp_regular=sp_regular.transpose("time", "freq", "dir")
    sp_regular['efth']=(('time', 'freq', 'dir'),sp_regular.efth.values)#/(1025*9.81))
    sp_regular['dir']=np.where(sp_regular.dir.values<0, sp_regular.dir.values+360, sp_regular.dir.values)
    sp_regular=sp_regular.sortby('dir')
    sp_regular
    m = np.reshape(np.sqrt(sp_regular.efth.values), (len(sp_regular.time), len(sp_regular.freq)*len(sp_regular.dir)))
    
    pca_fit = pk.load(open(os.path.join(path_data,'pca_fit.pkl'),'rb'))
    p=95
    
    APEV = np.cumsum(pca_fit.explained_variance_) / np.sum(pca_fit.explained_variance_)*100.0
    n_pcs=np.where(APEV>p)[0][0]
    print('n_pcs = ' + str(n_pcs))
    PCs_f = pca_fit.transform(m)
    
    #uw=sp_regular.Wspeed.values*np.cos(np.deg2rad(sp_regular.Wdir.values))
    #vw=sp_regular.Wspeed.values*np.sin(np.deg2rad(sp_regular.Wdir.values))
    uw=sp_regular.wind_x.values
    vw=sp_regular.wind_y.values
    matrix_rec= np.append(np.append(PCs_f[:,:n_pcs].T, [uw], axis=0), [vw], axis=0).T
    matrix_rec=np.where(np.isnan(matrix_rec)==True, 0, matrix_rec)
    
    #We load the mean and std values used to standarize the matrix in the fitting
    means=xr.open_dataset(os.path.join(path_data,'Sp_all_PCA_KMA.nc')).means.values
    stds=xr.open_dataset(os.path.join(path_data,'Sp_all_PCA_KMA.nc')).stds.values
    matrix_rec_norm=(matrix_rec-means)/stds
    
    np.shape(matrix_rec_norm) #Dimensions: Times x [n_pcs+2]
    
    kma_fit = pk.load(open(os.path.join(path_data,'kma_fit.pkl'),'rb'))
    predict_cluster = kma_fit.predict(matrix_rec_norm)
    
    
    kma_order=xr.open_dataset(os.path.join(path_data,'Sp_all_PCA_KMA.nc')).kma_order.values
    predict_cluster_sorted = np.zeros((len(predict_cluster),),)*np.nan
    for i in range(num_clusters):
        posc = np.where(predict_cluster == kma_order[i])
        predict_cluster_sorted[posc] = i
        
    #bmus=xr.open_dataset(os.path.join(path_data,'Sp_all_PCA_KMA.nc'))    
    return(sp_regular,predict_cluster_sorted)

def get_precomputed_simulations_and_create_nc_file(
        path_data,sp_regular,predict_cluster_sorted,nc_out_file
        ):
    swan_cluster = xr.open_dataset(os.path.join(path_data,'Adcirc_SWAN_Clusters.nc'))
    swan_time = sp_regular.time.values
    swan_x = swan_cluster.x.values
    swan_y = swan_cluster.y.values
    swan_hs = swan_cluster.Hs.values
    swan_tm01 = swan_cluster.Tm01.values
    swan_tp = swan_cluster.Tp.values
    swan_dir = swan_cluster.Dir.values
    swan_zeta = swan_cluster.zeta.values
    
    target_x = np.arange(np.min(swan_x)+0.15,np.max(swan_x)-0.15,0.0005)
    target_y = np.arange(np.min(swan_y)+0.1,np.max(swan_y)-0.1,0.0005)
    
    lgrid_x, lgrid_y = np.meshgrid(target_x,target_y)
    
    target_hs = np.ones((len(predict_cluster_sorted),len(target_y),len(target_x)))*-999
    target_tm01 = np.ones((len(predict_cluster_sorted),len(target_y),len(target_x)))*-999
    target_tp = np.ones((len(predict_cluster_sorted),len(target_y),len(target_x)))*-999
    target_dir = np.ones((len(predict_cluster_sorted),len(target_y),len(target_x)))*-999
    target_zeta = np.ones((len(predict_cluster_sorted),len(target_y),len(target_x)))*-999
    
    for ix in range(48,len(predict_cluster_sorted)):
        target_hs[ix,:,:] = sio.griddata(np.array([swan_x,swan_y]).T,swan_hs[np.int64(predict_cluster_sorted[ix]),:],(lgrid_x,lgrid_y))
        target_tm01[ix,:,:] = sio.griddata(np.array([swan_x,swan_y]).T,swan_tm01[np.int64(predict_cluster_sorted[ix]),:],(lgrid_x,lgrid_y))
        target_tp[ix,:,:] = sio.griddata(np.array([swan_x,swan_y]).T,swan_tp[np.int64(predict_cluster_sorted[ix]),:],(lgrid_x,lgrid_y))
        target_dir[ix,:,:] = sio.griddata(np.array([swan_x,swan_y]).T,swan_dir[np.int64(predict_cluster_sorted[ix]),:],(lgrid_x,lgrid_y))
        target_zeta[ix,:,:] = sio.griddata(np.array([swan_x,swan_y]).T,swan_zeta[np.int64(predict_cluster_sorted[ix]),:],(lgrid_x,lgrid_y))
        print(ix)
       
    df1 = xr.DataArray(target_hs, coords=[('time', swan_time), ('lat', target_y),('lon',target_x)])
    df1.name = 'Hs'
    df2 = xr.DataArray(target_tm01, coords=[('time', swan_time), ('lat', target_y),('lon',target_x)])
    df2.name = 'Tm'
    df3 = xr.DataArray(target_tp, coords=[('time', swan_time), ('lat', target_y),('lon',target_x)])
    df3.name = 'Tp'
    df4 = xr.DataArray(target_dir, coords=[('time', swan_time), ('lat', target_y),('lon',target_x)])
    df4.name = 'Dir'
    df5 = xr.DataArray(target_zeta, coords=[('time', swan_time), ('lat', target_y),('lon',target_x)])
    df5.name = 'Zeta'
    df = xr.merge([df1,df2,df3,df4,df5])
    df.to_netcdf(nc_out_file)
    return(swan_cluster)

def haversine_distance(lat1, lon1, lat2, lon2):
    p = 0.017453292519943295
    hav = 0.5 - cos((lat2-lat1)*p)/2 + cos(lat1*p)*cos(lat2*p) * (1-cos((lon2-lon1)*p)) / 2
    return 12742 * asin(sqrt(hav))

def extract_point_output(path_data,wb_lon,wb_lat,sp_regular,predict_cluster_sorted):
    swan_cluster = xr.open_dataset(os.path.join(path_data,'Adcirc_SWAN_Clusters.nc'))
    swan_x = swan_cluster.x.values
    swan_y = swan_cluster.y.values
    
    dis = np.zeros(np.size(swan_x))
    for ix in range(len(swan_x)):
        dis[ix] = haversine_distance(swan_y[ix],swan_x[ix],wb_lat,wb_lon)
    posix = np.where(dis == np.min(dis))
    posix = posix[0][0]
    
    swan_time = sp_regular.time.values
    swan_x_pt = swan_cluster.x.values[posix]
    swan_y_pt = swan_cluster.y.values[posix]
    
    swan_hs_pt = np.zeros(len(swan_time))
    swan_tm01_pt = np.zeros(len(swan_time))
    swan_tp_pt = np.zeros(len(swan_time))
    swan_dir_pt = np.zeros(len(swan_time))
    swan_zeta_pt = np.zeros(len(swan_time))
    for ix in range(len(swan_time)):
        swan_hs_pt[ix] = swan_cluster.Hs.values[np.int64(predict_cluster_sorted[ix]),posix]
        swan_tm01_pt[ix] = swan_cluster.Tm01.values[np.int64(predict_cluster_sorted[ix]),posix]
        swan_tp_pt[ix] = swan_cluster.Tp.values[np.int64(predict_cluster_sorted[ix]),posix]
        swan_dir_pt[ix] = swan_cluster.Dir.values[np.int64(predict_cluster_sorted[ix]),posix]
        swan_zeta_pt[ix] = swan_cluster.zeta.values[np.int64(predict_cluster_sorted[ix]),posix]
    return(swan_time,swan_x_pt,swan_y_pt,swan_hs_pt,swan_tm01_pt,swan_tp_pt,swan_dir_pt,swan_zeta_pt)

def make_ts_fig(time_waves,twl_nearshore,tide,time_tide_minute,tide_minute,mx_lvl1,mx_lvl2,outfile_name):
    time_waves_pd = pd.to_datetime(time_waves)
    
    xor = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_waves]
    xnew = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_tide_minute]
    f = sio.interp1d(xor,twl_nearshore,kind='cubic',fill_value="extrapolate")
    twl_min=f(xnew)
    
    limits = [-1,np.max(mx_lvl2)+0.25]
    dayFmt = mdates.DateFormatter('%d-%m-%Y %H:%M')
    
    #-- differentiate to calculate high and low tides
    diff = np.zeros_like(time_tide_minute, dtype=np.float64)
    #-- forward differentiation for starting point
    diff[0] = tide_minute[1] - tide_minute[0]
    #-- backward differentiation for end point
    diff[-1] = tide_minute[-1] - tide_minute[-2]
    #-- centered differentiation for all others
    diff[1:-1] = (tide_minute[2:] - tide_minute[0:-2])/2.0
    htindex, = np.nonzero((np.sign(diff[0:-1]) >= 0) & (np.sign(diff[1:]) < 0))
    # ltindex, = np.nonzero((np.sign(diff[0:-1]) <= 0) & (np.sign(diff[1:]) > 0))

    #plt.figure(figsize=(10,6))
    f, ax = plt.subplots(figsize=(12,6))
    
    #ax.axhspan(mx_lvl1,mx_lvl2,facecolor='yellow',alpha=0.5)
    #ax.axhspan(mx_lvl2,8,color='r',alpha=0.5)
    
    ax.plot(time_waves_pd[48:],np.matlib.repmat(mx_lvl2,len(time_waves_pd[48:]),1),color='firebrick')
    ax.text(time_waves_pd[-55], mx_lvl2+0.05, 'Moderate flood threshold',color='firebrick',size=11)
    ax.plot(time_waves_pd[48:],np.matlib.repmat(mx_lvl1,len(time_waves_pd[48:]),1),color='darkorange')
    plt.text(time_waves_pd[-50], mx_lvl1+0.05, 'Minor flood threshold',color='darkorange',size=11)
    
    ax.fill_between(time_waves_pd[48:],twl_nearshore[48:],-2,color='lightsteelblue',label='Nearshore TWL')
    ax.fill_between(time_waves_pd[48:],tide[48:],-2,color='royalblue',label='Offshore Tide')
    for h in range(4,len(htindex)-1):
            text=time_tide_minute[htindex[h]].strftime("%H:%M")
            plt.plot(time_tide_minute[htindex[h]],twl_min[htindex[h]],"v",color= 'gray',markersize=5)
            plt.text(time_tide_minute[htindex[h]]-dt.timedelta(hours=3),twl_min[htindex[h]]+0.08,text,color= 'gray',size=15) 
    ax.set_xlim(time_waves_pd[48],time_waves_pd[-1:])
    ax.set_ylim(limits)
    ax.set_ylabel('Water level (m)')
    ax.grid(True,alpha = 0.7)
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour = [6,12,18]))
    ax.xaxis.set_major_formatter(dayFmt)
    f.autofmt_xdate()
    ax.legend()
    f.savefig(outfile_name)
    plt.close(f)
    return()


def readvars_nc(nc_fname,varname):
    ncc = nc.Dataset(nc_fname,'r') 
    if varname == 'time':
        time_var = ncc.variables[varname] 
        time_or = nc.num2date(ncc[varname],time_var.units, time_var.calendar)
        time_str=[time_or[i].strftime('%Y-%m-%d %H:%M') for i in range(0,len(time_or))]
        var = [dt.datetime.strptime(time_str[i],'%Y-%m-%d %H:%M') for i in range(0,len(time_str))]
    else:
        var = np.array(ncc[varname])
    ncc.close()
    return(var)

def make_risk_csv_file(nearshore_twl,mx_lvl1,mx_lvl2,lon,lat,risk_outfile_name):
    output_twl_nearshore = nearshore_twl[48:]
    if output_twl_nearshore.max() >= mx_lvl2:
        risk_cat = 2
    elif output_twl_nearshore.max() >= mx_lvl1:
        risk_cat = 1
    else:
        risk_cat = 0
    
    risk_output = np.array((lon,lat,risk_cat))
    risk_output = np.reshape(risk_output,(len(risk_output),1)).T
    np.savetxt(risk_outfile_name, risk_output, delimiter=",")
    return()

def make_lagoon_forecast(now,path):
    path_data=os.path.join(path, 'Data')
    location_path=os.path.join(path, 'Funafuti_Lagoonside_points/Funafuti_LagoonSide_Points_v3.csv')    
    transects = pd.read_csv(location_path)

    fl_name_tides = os.path.join(path, '../tmp/tide_hourly_Fongafale.nc')
    fl_name_tides_min = os.path.join(path, '../tmp/tide_minute_Fongafale.nc')
    fl_name_sla = os.path.join(path, '../tmp/sla_hourly_Fongafale.nc')
    time_sla = readvars_nc(fl_name_sla,'time')
    sla = readvars_nc(fl_name_sla,'SLA')/100
    time_tide = readvars_nc(fl_name_tides,'time')
    tide = readvars_nc(fl_name_tides,'tide')/100
    time_tide_minute = readvars_nc(fl_name_tides_min,'time')
    tide_minute = readvars_nc(fl_name_tides_min,'tide')/100
    
    create_superpoint_nc(now)    
    sp_regular,predict_cluster_sorted=determine_clusters(now,path_data)
    
    #nc_out_file = os.path.join(path, '../runs/'+now.strftime('%Y%m%d%H')+'/results/test.nc')
    #swan_cluster = get_precomputed_simulations_and_create_nc_file(path_data,sp_regular,predict_cluster_sorted,nc_out_file)
    
    xor = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_sla]
    xnew = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_tide]
    f = sio.interp1d(xor,sla,kind='cubic',fill_value="extrapolate")
    sla_1d_new=f(xnew)
    Pt_SLA_Tide = tide + sla_1d_new
    print('Creating lagoon output')
    for i in range(len(transects)):
        thresholds_path = os.path.join(path, 'Thresholds/t_' + str(i) + '.csv')
        thresholds = np.loadtxt(thresholds_path,delimiter=',',skiprows=0)
        mx_lvl1 = thresholds[0,1]
        mx_lvl2 = thresholds[72,1]
        # mx_lvl1 = 3
        # mx_lvl2 = 4
        outfile_name = '../inundation/Figures/Funafutilagoon_t_' + str(i) + '_forecast.png'
        risk_outfile_name = '../inundation/Flood_risk/Funafutilagoon_t_' + str(i) + '_risk_level.csv'

        t_lat = transects.X[i]
        t_lon = transects.Y[i]#lat and lon are flipped in the point csv file
        swan_time,swan_x_pt,swan_y_pt,swan_hs_pt,swan_tm01_pt,swan_tp_pt,swan_dir_pt,swan_zeta_pt = extract_point_output(path_data,t_lon,t_lat,sp_regular,predict_cluster_sorted)
        ix = np.isnan(swan_zeta_pt)
        swan_zeta_pt[ix] = 0
        nearshore_twl = Pt_SLA_Tide + swan_zeta_pt
        
        time_waves = [pd.Timestamp(swan_time[j]).to_pydatetime() for j in range(0,len(swan_time))]
        make_ts_fig(time_waves,nearshore_twl,tide,time_tide_minute,tide_minute,mx_lvl1,mx_lvl2,outfile_name)
        make_risk_csv_file(nearshore_twl,mx_lvl1,mx_lvl2,t_lon,t_lat,risk_outfile_name)
        print(str(i+1) + ' / ' + str(len(transects)))
    return()

#import sys

#path=r'F:\Adcirc_SWAN\Forecast_system\Funafuti_lagoon'
#sys.path.insert(1, path)
#now = dt.datetime(2022,1,16,6)

#make_lagoon_forecast(now,lagoon_path)

# efth_cluster = xr.open_dataset(os.path.join(path_data,'Centroids_info_2000.nc'))

# from matplotlib import pyplot as plt
# cm = 'jet'
# plt.figure()
# plt.subplot(1,2,1)
# plt.pcolor(sp_regular.dir.values,sp_regular.freq.values,
#             sp_regular.efth.values[180,:,:],
#             cmap=cm, vmin=0, vmax=0.1)
# plt.colorbar()
# plt.subplot(1,2,2)
# plt.pcolor(
#     efth_cluster.dir.values,efth_cluster.freq.values,
#     efth_cluster.efth.values[np.int64(predict_cluster_sorted[180]),:,:],
#     cmap=cm, vmin=0, vmax=0.1)
# plt.colorbar()





    