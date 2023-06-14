# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 13:59:09 2023

@author: moritzw
"""

import os, glob
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import subprocess
from datetime import timedelta


def convert_grib_2_nc(grb_fl_name,nc_fl_name):
# function to convert grib2 files to netCDF, it uses the java script "toolsUI-5.4.1.jar" that needs to be placed in the same folder where this script is being called    
    subprocess.call(["java", "-Xmx512m", "-cp", "toolsUI-5.4.1.jar", "ucar.nc2.dataset.NetcdfDataset", "-in", grb_fl_name, "-out", nc_fl_name], shell=True,)
    print(grb_fl_name + ' converted to ' + nc_fl_name)
    return()


def download_all_regional_data_grb_and_convert_2_nc(mydate,Tcycle,dt,end_tt,leftlon,rightlon,toplat,bottomlat,grb_out,nc_out):

# function to download the winds
# mydate= current date
# Tcycle= forecast run, 00,06,12,18
# dt= time step, usually 3 hours
# end_dt= end of the forecast period
# leftlon, rightlon,toplat, bottomlat= limits of the area to be downloaded
# grb_out= root name of grib2 files
# nc_out= root name of netcdf file
    
    # define parameters for the connection with the server
    session = requests.Session()
    retry = Retry(connect=50, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
      
    # there are independent files for each time step, we need to go through all the times in the forecast period
    for i in range(0,end_tt,dt):
        
        # url to define the time and region to download
        url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?' + \
        'file=gfs.t' + Tcycle + 'z' + \
        '.pgrb2.0p25.f' + "{0:0>3}".format(i) + \
        '&lev_10_m_above_ground=on&lev_surface=on&lev_mean_sea_level=on&var_PRMSL=on&var_UGRD=on&var_VGRD=on&var_VIS=on&var_GUST=on&var_TMP=on' + \
        '&subregion=&leftlon=' + leftlon + \
        '&rightlon=' + rightlon + \
        '&toplat=' + toplat + \
        '&bottomlat=' + bottomlat + \
        '&dir=%2Fgfs.' + mydate + '%2F' + Tcycle + '%2Fatmos'
        
        session.mount(url, adapter)  
        r=session.get(url)  
        

        grb_outX = grb_out + '_' + "{0:0>3}".format(i) +'.grib2'
        nc_outX = nc_out + '_' + "{0:0>3}".format(i) +'.nc'

        # write grib2 on disk
        open(grb_outX, 'wb').write(r.content)
        print('Grib file downloaded and stored as ' + grb_outX) 
	# convert to netcdf       
        convert_grib_2_nc(grb_outX,nc_outX)
    session.close()    
    return()


def download_NCEP(now):
    # now= date and hour of the NOAA/NCEP run to download, it is a datetime object but this can be changed to our needs. I define this date in a previous script that checks the latest available run on the server
    Tcycle = str(now.hour).zfill(2)
    # time interval of the forecast, usually 3 hours
    wave_dt = 3
    wind_dt = 3

    # from the current time downlad the next 180 (+1) hours (7.5 days)
    time_length = 181
    # I also dowload the first 21 hours of the previous day to spin up the model (1 day), if you use hotstart you can negret this step
    htime_length = 22
    # coordinates of the region to be downloades
    leftlon = '140.0'
    rightlon = '210.0'
    toplat = '30.0'
    bottomlat = '-60.0'
      
    ###############################################################################
    # mydate is a string 20210902
    mydate = now.strftime("%Y") + now.strftime("%m") + now.strftime("%d")

    wnd_grb_out = '../Regional_tmp/out'
    wnd_nc_out = '../Regional_tmp/out'
    # download wind and waves for the forecast period (7.5 days)
    download_all_regional_data_grb_and_convert_2_nc(mydate,Tcycle,wind_dt,time_length,leftlon,rightlon,toplat,bottomlat,wnd_grb_out,wnd_nc_out)

    # remove all the grib2files
    for filename in glob.glob("../tmp/*.grib2"):
        os.remove(filename) 
    
    