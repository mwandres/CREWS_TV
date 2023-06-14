# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 09:40:26 2021

@author: 
"""
import os
import netCDF4
import numpy as np
from netCDF4 import num2date
import pandas as pd
import create_spec_from_partitions

# Functions
def read_netcdf(nc_fname,dt,length_dt):
# Function to read wave partitions contained in the netCDF files from step1_download_NCEP.py, gives back wave partitions for the whole tile. 
# This function also corrects wave height partitions to make their summation to match with the total Hs produced by WWIII before partitioning the spectra

    # read wave partition for hour=0
    nc_fnameX = nc_fname + '_' + "{0:0>3}".format(0) +'.nc'
    nc = netCDF4.Dataset(nc_fnameX)
    # kk = np.array(nc['ordered_sequence_of_data'])
    ln = np.array(nc['lon'])
    lt = np.array(nc['lat'])
    tt = np.array(nc['time'])
    Hs_ST = np.nan_to_num(np.array(nc['Significant_height_of_combined_wind_waves_and_swell_surface'])[:,:,:])
    # Tm_ST = np.nan_to_num(np.array(nc['Primary_wave_mean_period_surface'])[:,:,:])
    # Dir_ST = np.nan_to_num(np.array(nc['Primary_wave_direction_surface'])[:,:,:])
    Hs_S1 = np.nan_to_num(np.array(nc['Significant_height_of_swell_waves_ordered_sequence_of_data'])[:,0,:,:])
    Hs_S2 = np.nan_to_num(np.array(nc['Significant_height_of_swell_waves_ordered_sequence_of_data'])[:,1,:,:])
    Hs_S3 = np.nan_to_num(np.array(nc['Significant_height_of_swell_waves_ordered_sequence_of_data'])[:,2,:,:])
    Hs_W = np.nan_to_num(np.array(nc['Significant_height_of_wind_waves_surface'])[:,:,:])   
    Tm_S1 = np.nan_to_num(np.array(nc['Mean_period_of_swell_waves_ordered_sequence_of_data'])[:,0,:,:])
    Tm_S2 = np.nan_to_num(np.array(nc['Mean_period_of_swell_waves_ordered_sequence_of_data'])[:,1,:,:])
    Tm_S3 = np.nan_to_num(np.array(nc['Mean_period_of_swell_waves_ordered_sequence_of_data'])[:,2,:,:])
    Tm_W = np.nan_to_num(np.array(nc['Mean_period_of_wind_waves_surface'])[:,:,:])
    Dir_S1 = np.nan_to_num(np.array(nc['Direction_of_swell_waves_ordered_sequence_of_data'])[:,0,:,:])
    Dir_S2 = np.nan_to_num(np.array(nc['Direction_of_swell_waves_ordered_sequence_of_data'])[:,1,:,:])
    Dir_S3 = np.nan_to_num(np.array(nc['Direction_of_swell_waves_ordered_sequence_of_data'])[:,2,:,:])
    Dir_W = np.nan_to_num(np.array(nc['Direction_of_wind_waves_surface'])[:,:,:])    
   # read wave partition from hour=3 to hour=180 every 3 hours
    for i in range(dt,length_dt,dt):
        nc_fnameX = nc_fname + '_' +  "{0:0>3}".format(i) +'.nc'
        nc = netCDF4.Dataset(nc_fnameX)
        tt = np.append(tt,np.array(nc['time']),axis=0)
        Hs_ST = np.nan_to_num(np.append(Hs_ST,np.array(nc['Significant_height_of_combined_wind_waves_and_swell_surface'])[:,:,:],axis=0))
        Hs_S1 = np.nan_to_num(np.append(Hs_S1,np.array(nc['Significant_height_of_swell_waves_ordered_sequence_of_data'])[:,0,:,:],axis=0))
        Hs_S2 = np.nan_to_num(np.append(Hs_S2,np.array(nc['Significant_height_of_swell_waves_ordered_sequence_of_data'])[:,1,:,:],axis=0))
        Hs_S3 = np.nan_to_num(np.append(Hs_S3,np.array(nc['Significant_height_of_swell_waves_ordered_sequence_of_data'])[:,2,:,:],axis=0))
        Hs_W = np.nan_to_num(np.append(Hs_W,np.array(nc['Significant_height_of_wind_waves_surface'])[:,:,:],axis=0))
        Tm_S1 = np.nan_to_num(np.append(Tm_S1,np.array(nc['Mean_period_of_swell_waves_ordered_sequence_of_data'])[:,0,:,:],axis=0))
        Tm_S2 = np.nan_to_num(np.append(Tm_S2,np.array(nc['Mean_period_of_swell_waves_ordered_sequence_of_data'])[:,1,:,:],axis=0))
        Tm_S3 = np.nan_to_num(np.append(Tm_S3,np.array(nc['Mean_period_of_swell_waves_ordered_sequence_of_data'])[:,2,:,:],axis=0))
        Tm_W = np.nan_to_num(np.append(Tm_W,np.array(nc['Mean_period_of_wind_waves_surface'])[:,:,:],axis=0))
        Dir_S1 = np.nan_to_num(np.append(Dir_S1,np.array(nc['Direction_of_swell_waves_ordered_sequence_of_data'])[:,0,:,:],axis=0))
        Dir_S2 = np.nan_to_num(np.append(Dir_S2,np.array(nc['Direction_of_swell_waves_ordered_sequence_of_data'])[:,1,:,:],axis=0))
        Dir_S3 = np.nan_to_num(np.append(Dir_S3,np.array(nc['Direction_of_swell_waves_ordered_sequence_of_data'])[:,2,:,:],axis=0))
        Dir_W = np.nan_to_num(np.append(Dir_W,np.array(nc['Direction_of_wind_waves_surface'])[:,:,:],axis=0))
        print(nc_fnameX)
        
    # determine calbration coeficients
    factor= Hs_ST/np.sqrt(Hs_S1**2+Hs_S2**2+Hs_S3**2+Hs_W**2)
    # apply calbration coeficients
    Hs_S1=Hs_S1*factor
    Hs_S2=Hs_S2*factor
    Hs_S3=Hs_S3*factor
    Hs_W=Hs_W*factor
    
    
    dates = num2date(tt[:], units=nc['time'].units, calendar = nc['time'].calendar)
    return(ln,lt,dates,Hs_S1,Hs_S2,Hs_S3,Hs_W,Tm_S1,Tm_S2,Tm_S3,Tm_W,Dir_S1,Dir_S2,Dir_S3,Dir_W)

def makeSWANspec(out_name,time,fff,ddd,EE,x,y,spherical):
    print('Writing spectral file ' + out_name)
    #    function makeSWANspec(filename,time,fff,ddd,EE,x,y,spherical)
    #%create a timeseries 2d spectral file for SWAN
    #%makeSWANspec(filename,time,fff,ddd,EE,x,y,spherical)
    #%
    #%filename: filename for output spectral file
    #%time:time vector for the timeseries
    #%frequency vector for the 2d spectra (in Hz)
    #%directional vector in Deg
    #%EE var dens spectra time series as EE(ddd,fff,time)
    #%x,y are the location of the timeseries
    #%spherical is a switch to say whether x and y are lon and lat (=1) or
    #%catesian coordiantes (=0)
    #%CB 09-09-2013
    output =open(out_name,'w')
    s='SWAN  1\n';
    output.write(s);

    if len(time)>1:
        s='TIME\n\t1\n';
        output.write(s);
        
    if spherical>0:
        s='LONLAT\n';
        output.write(s);
    else:
        s='LOCATIONS\n';
        output.write(s);
    
    s='\t1\n%f\t%f\n' %(x,y);
    output.write(s);
    s='AFREQ\n%d\n' %len(fff);
    output.write(s);
    
    for i in range(len(fff)):
        s='%f\n' %fff[i]
        output.write(s)
    
    s='NDIR\n%d\n' %len(ddd);
    output.write(s);

    for i in range(len(ddd)):
        s='%f\n' %ddd[i];
        output.write(s);

    s='QUANT\n';
    output.write(s);
    s='\t1\n';
    output.write(s);    
    s='VaDens\n';
    output.write(s);
    s='m2/Hz/degr\n';
    output.write(s);
    s='\t-0.9900E+02\n';
    output.write(s);

    for n in range(len(time)):
        if len(time)>1:
            s='%s\n' %time[n].strftime('%Y%m%d.%H%M%S');#'yyyymmdd.HHMMSS'));
            output.write(s);
                        
        s='FACTOR\n';
        output.write(s);        
        s='0.000001\n';
        output.write(s);
        
        for fi in range(len(fff)):
            for di in range(len(ddd)):
                s='%d\t' %round(EE[n,fi,di]/0.000001);
                output.write(s);

            s='\n';
            output.write(s);
    output.close()
    return()

def extract_entire_time_series(time,x,y,freqs,theta,Hs_S1,Hs_S2,Hs_S3,Hs_W,Tm_S1,Tm_S2,Tm_S3,Tm_W,Dir_S1,Dir_S2,Dir_S3,Dir_W):
# Function to extract time series of wave partitions at the required location,
# determine directional spreading accoring to a relation obtained in the area
# returns the 2D spectra for all times at one location
# x,y are the indices of the point
# freqs frequency vector for the 2d spectra (in Hz)
# theta, directional vector in Deg
# Hs_S... wave partition parameters

    MAT = np.zeros((len(time),len(freqs),len(theta)))
    for ttt in range(len(time)):
      
        Hs1 = Hs_S1[ttt,x,y]
        Dir1 = Dir_S1[ttt,x,y]
        Tp1 = Tm_S1[ttt,x,y]
        Hs2 = Hs_S2[ttt,x,y]
        Dir2 = Dir_S2[ttt,x,y]
        Tp2 = Tm_S2[ttt,x,y]
        Hs3 = Hs_S3[ttt,x,y]
        Dir3 = Dir_S3[ttt,x,y]
        Tp3 = Tm_S3[ttt,x,y]
        Hs4 = Hs_W[ttt,x,y]
        Dir4= Dir_W[ttt,x,y]
        Tp4 = Tm_W[ttt,x,y]
        
        phs = np.squeeze([Hs1,Hs2,Hs3,Hs4]).T   
        # hstotal=np.sqrt(np.sum(phs.T**2))
        ptp = np.squeeze([Tp1,Tp2,Tp3,Tp4]).T
        pth = np.squeeze([Dir1,Dir2,Dir3,Dir4]).T
        

        #Best Tp-psi fit for directional spreading (in degrees) in the area obtained from CAWCR reanalysis, directional spreading values shpuld be defined otherwise
        #psi=np.array([20,20,20,20])
        psi=(7.993e-05*ptp**5)-(0.007424*ptp**4)+(0.243*ptp**3)-(3.416*ptp**2)+(18.01*ptp)-0.003026
        psi[ptp>20]=5.7
        psi[ptp<2]=30
        # define peak enhancement factor of the JONSWAP spectra
        pgamma = 3.3
        # determine wave spectra for all times
        Energy, freqsMatrix, dirsMatrix = create_spec_from_partitions.freqdir_spectrum(freqs, theta, phs, ptp, pgamma, psi, pth)
        MAT_combined = np.squeeze(Energy.T)
        MAT[ttt,:,:] = MAT_combined

    return(MAT)

###############################################################################
def make_waves(now):
    
    print('Generating spectral wave forcing')  
    
    # define the output folder were SWAN will be run
    out_name = '../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H") 

    ###############################################################################
    #freqs = np.arange(0.03,1.,.01)
    # logarithmic scale definition inf frequencise as it is in SWAN manual,  increases resolution on the long periods, larger bins on the short waves 
    freqs = 0.0373*np.power(1.1,np.arange(0,32,1))
    theta = np.arange(0,360,10)
    ###############################################################################
   
    # directory of the nc files
    wave_nc = '../tmp/wave_tmp'
    h_wave_nc = '../tmp/h_wave_tmp'
    hh_wave_nc = '../tmp/hh_wave_tmp'
     
    wave_dt = 3
    time_length = 181
    h_time_length = 22
     
    # read wave nc files
    ln,lt,tt,Hs_S1,Hs_S2,Hs_S3,Hs_W,Tm_S1,Tm_S2,Tm_S3,Tm_W,Dir_S1,Dir_S2,Dir_S3,Dir_W = read_netcdf(wave_nc,wave_dt,time_length) 
    ln,lt,h_tt,h_Hs_S1,h_Hs_S2,h_Hs_S3,h_Hs_W,h_Tm_S1,h_Tm_S2,h_Tm_S3,h_Tm_W,h_Dir_S1,h_Dir_S2,h_Dir_S3,h_Dir_W = read_netcdf(h_wave_nc,wave_dt,h_time_length) 
    ln,lt,hh_tt,hh_Hs_S1,hh_Hs_S2,hh_Hs_S3,hh_Hs_W,hh_Tm_S1,hh_Tm_S2,hh_Tm_S3,hh_Tm_W,hh_Dir_S1,hh_Dir_S2,hh_Dir_S3,hh_Dir_W = read_netcdf(hh_wave_nc,wave_dt,h_time_length) 
     
    #%%
      
    try:
        os.mkdir(out_name)
    except OSError as error:
        print(error)    
  
     # read the location of the boundary points of the SWAN unestructured mesh
    f14bnds = pd.read_csv("../extras/f14_boundary_points.csv",sep=",")
    f14_xb = np.array(f14bnds["x_f14_bnd"])
    f14_yb = np.array(f14bnds["y_f14_bnd"])
    f14_id = np.array(f14bnds["adcirc_index"])
  
    
    # find the indices of the closest point in the wave grid to each boundary point of SWAN mesh
    yv, xv = np.meshgrid(lt,ln, sparse=False, indexing='ij')   
    lon=np.reshape(xv, (np.size(xv), 1)) 
    lat=np.reshape(yv, (np.size(yv), 1))
     
    Wf_index=[]
    for k in range(len(f14_xb)):
        dist=np.squeeze(np.sqrt((lon-f14_xb[k])**2+(lat-f14_yb[k])**2))
        kk=np.argsort(dist)[0]
        Wf_index.append(kk) 
       
     # this part smooths the wave spectra along the boundary, needed when you work with unstructured SWAN meshes where boundary points are  defined independently  
    shape=np.shape(xv)
    Wf_index=np.asarray( Wf_index)
     
    for l in range(len(f14_xb)):
        
        pos=np.arange(l-1,l+2)
        pos[pos<0]=pos[pos<0]+len(f14_xb)
        pos[pos>len(f14_xb)-1]=pos[pos>len(f14_xb)-1]-len(f14_xb)
        
        out_nameXX=out_name+'/Pto_' + str(f14_id[l]) + '.sp2'   
        [i, j] = np.unravel_index(Wf_index[pos], shape) 
         
        indices=np.vstack([i,j]).T
        uind=np.unique(indices,axis=0)
        
        MATT_hh=np.zeros((np.size(hh_tt),np.size(freqs),np.size(theta)))
        MATT_h=np.zeros((np.size(h_tt),np.size(freqs),np.size(theta)))
        MATT_f=np.zeros((np.size(tt),np.size(freqs),np.size(theta)))
        
  
        
        for m in range(len(uind)):
            
            iii=uind[m,0]
            jjj=uind[m,1]
           
            MAT_hh = extract_entire_time_series(hh_tt,iii,jjj,freqs,theta,hh_Hs_S1,hh_Hs_S2,hh_Hs_S3,hh_Hs_W,hh_Tm_S1,hh_Tm_S2,hh_Tm_S3,hh_Tm_W,hh_Dir_S1,hh_Dir_S2,hh_Dir_S3,hh_Dir_W)
            MATT_hh=MATT_hh+MAT_hh
            MAT_h = extract_entire_time_series(h_tt,iii,jjj,freqs,theta,h_Hs_S1,h_Hs_S2,h_Hs_S3,h_Hs_W,h_Tm_S1,h_Tm_S2,h_Tm_S3,h_Tm_W,h_Dir_S1,h_Dir_S2,h_Dir_S3,h_Dir_W)
            MATT_h=MATT_h+MAT_h
            MAT_f = extract_entire_time_series(tt,iii,jjj,freqs,theta,Hs_S1,Hs_S2,Hs_S3,Hs_W,Tm_S1,Tm_S2,Tm_S3,Tm_W,Dir_S1,Dir_S2,Dir_S3,Dir_W)
            MATT_f=MATT_f+MAT_f
            
        TT=np.concatenate((hh_tt,h_tt,tt), axis=0) 
        MATT=np.concatenate((MATT_hh,MATT_h,MATT_f), axis=0)/len(uind)
        
        # create SWAN boundary file for each point in the boundary, this can be used to generate a file with all the points for regular SWAN grids
        makeSWANspec(out_nameXX,TT,freqs,theta,MATT,f14_xb[l],f14_yb[l],1)
        # time_str=[TT[i].strftime('%Y-%m-%d %H:%M') for i in range(0,len(tt))]
    print('Wave forcing generated in: ' + out_name)