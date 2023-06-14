"""
Created on Tue Dec 14 20:54:13 2021

@author: moritzw
"""
import scipy.io as sio
import numpy as np
import numpy.matlib
import datetime as dt
import geopy.distance
import netCDF4 as nc
import pandas as pd
from rbf import RBF_Interpolation, Normalize
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import datetime
from scipy.interpolate import interp1d

def load_SWAN_data(fl_name,var_name):
     if var_name == 'time':
        time= True
        avname = 'Hs'
     else :
        time = False
        avname = var_name
     infos = sio.whosmat(fl_name)
     names, sizes, types  = zip(*infos)
     vnames = []
     for n in names:
        if avname in n[0:len(avname):1]:
            vnames.append(n)
     vnames.sort()
     if time :
        data = []
        for i in range(len(vnames)): 
            aux = vnames[i][len(vnames[0])-15:]
            data.append(aux)
     else:    
        MAT = sio.loadmat(fl_name)         
        data = np.zeros([len(vnames),len(MAT[vnames[0]].T)])
        for i in range(len(vnames)):            
            data[i,:] = MAT[vnames[i]]   
     return(data)
 
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
 
def find_coordinate_ix(target_lon,target_lat,lon_array,lat_array):
    coords_1 = np.array((target_lat, target_lon))
    dist = geopy.distance.distance(coords_1, np.array((lat_array[0],lon_array[0]))).km
    for i in range(1,len(lon_array)):
        dist = np.append(dist,geopy.distance.distance(coords_1, np.array((lat_array[i],lon_array[i]))).km)
    idx = dist.argmin()
    return(idx,lon_array[idx],lat_array[idx])   

def make_ts_fig(time_waves,twl_nearshore,tide,time_tide_minute,tide_minute,mx_lvl1,mx_lvl2,outfile_name):
    time_waves_pd = pd.to_datetime(time_waves)
    
    xor = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_waves]
    xnew = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_tide_minute]
    f = interp1d(xor,twl_nearshore,kind='cubic',fill_value="extrapolate")
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

def make_risk_csv_file(output_twl_nearshore,mx_lvl1,mx_lvl2,lon,lat,risk_outfile_name):
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

def inundation_forecast(now):
    fl_name_waves = '../runs/'+now.strftime('%Y%m%d%H')+'/output.mat'
    fl_name_tides = '../tmp/tide_hourly_Nui.nc'
    fl_name_tides_min = '../tmp/tide_minute_Nui.nc'
    fl_name_sla = '../tmp/sla_hourly_Nui.nc'
    
    Time = load_SWAN_data(fl_name_waves,"time")
    time_waves = [dt.datetime.strptime(Time[i],'%Y%m%d_%H%M%S') for i in range(0,len(Time))]
    Xp = load_SWAN_data(fl_name_waves,"Xp")
    Yp = load_SWAN_data(fl_name_waves,"Yp")
    Hs = load_SWAN_data(fl_name_waves,"Hsig")
    Tp = load_SWAN_data(fl_name_waves,"RTpeak")
    Tm = load_SWAN_data(fl_name_waves,"Tm02")
    Dir = load_SWAN_data(fl_name_waves,"Dir_")
    Windv_x = load_SWAN_data(fl_name_waves,"Windv_x")
    Windv_y = load_SWAN_data(fl_name_waves,"Windv_y")
    time_sla = readvars_nc(fl_name_sla,'time')
    sla = readvars_nc(fl_name_sla,'SLA')/100
    time_tide = readvars_nc(fl_name_tides,'time')
    tide = readvars_nc(fl_name_tides,'tide')/100
    time_tide_minute = readvars_nc(fl_name_tides_min,'time')
    tide_minute = readvars_nc(fl_name_tides_min,'tide')/100
    
    
    transects_csv = '../inundation/RBF/Nui_forecast/Nui/Nui_Oceanside_Profilesv4.csv'
    transects = pd.read_csv(transects_csv,delimiter=',',skiprows=0)
    transect_labels =pd.unique(transects.LINE_ID)
    transect_labels = transect_labels.astype(int)
    
    for tt in range(len(transect_labels)):
        fig_outfile_name = '../inundation/Figures/Nui_t_' + str(transect_labels[tt]) + '_forecast.png'
        risk_outfile_name = '../inundation/Flood_risk/Nui_t_' + str(transect_labels[tt]) + '_risk_level.csv'
        
        coordinate_fl = '../inundation/RBF/Nui_forecast/t_' + str(transect_labels[tt]) + '_coordinates.csv'
        training_data_path = '../inundation/RBF/Nui_forecast/t_' + str(transect_labels[tt]) + '_results.csv'
        rbf_coeff_path = '../inundation/RBF/Nui_forecast/t_' + str(transect_labels[tt]) + '_rbf_coefficients_max_TWL_and_inun_ext.pkl'
        thresholds_path = '../inundation/RBF/Nui_forecast/t_' + str(transect_labels[tt]) + '_thresholds.csv'
        
        coords = np.loadtxt(coordinate_fl,delimiter=',',skiprows=0)
        target_x = coords[0]
        target_y = coords[1]
        
        idx,lon,lat = find_coordinate_ix(target_x, target_y, Xp[0,:], Yp[0,:])
        Pt_Hs = Hs[:,idx]
        Pt_Tp = Tp[:,idx]
        Pt_Dir = Dir[:,idx]
        #sla_1d = np.squeeze(sla)
        xor = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_sla]
        xnew = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_waves]
        f = interp1d(xor,sla,kind='cubic',fill_value="extrapolate")
        sla_1d_new=f(xnew)
        Pt_SLA_Tide = tide + sla_1d_new
        Pt_target_data = np.array((Pt_Hs,Pt_Tp,Pt_Dir,Pt_SLA_Tide)).T
        
        training_data_csv = np.loadtxt(training_data_path,delimiter=',',skiprows=0)
        training_dataset = training_data_csv[:,0:4]
        # subset - scalar / directional indexes
        ix_scalar = [0,1,3]      # scalar (hs, tp)
        ix_directional = [2]   # directional (dir)
        ##
        df = pd.read_pickle(rbf_coeff_path)
        
        # normalize subset and dataset
        subset_norm, mins, maxs = Normalize(
            training_dataset, ix_scalar, ix_directional)
        
        dataset_norm, _, _ = Normalize(
            Pt_target_data, ix_scalar, ix_directional,mins,maxs)
        
        opt_sigma = df.opt_sigma[0]
        rbf_coeff = df.rbf_coeff[0]
        output_twl_nearshore = RBF_Interpolation(opt_sigma, rbf_coeff, subset_norm.T, dataset_norm.T)
        output_twl_nearshore = np.nan_to_num(output_twl_nearshore)
        
        thresholds = np.loadtxt(thresholds_path,delimiter=',',skiprows=0)
        mx_lvl1 = thresholds[0,1]
        mx_lvl2 = thresholds[72,1]
        
        make_ts_fig(time_waves,output_twl_nearshore,tide,time_tide_minute,tide_minute,mx_lvl1,mx_lvl2,fig_outfile_name)
        make_risk_csv_file(output_twl_nearshore,mx_lvl1,mx_lvl2,coords[2],coords[3],risk_outfile_name)
    return()
    
# now = dt.datetime(2021,12,22,6)
# inundation_forecast(now)