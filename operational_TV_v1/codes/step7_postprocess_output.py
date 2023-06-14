# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 08:51:33 2021

@author: judithg
"""
# import math
from adcircpy import AdcircMesh
import pandas as pd
import numpy as np
import datetime as dt
from matplotlib import pyplot as plt
import matplotlib.tri as tri
import scipy.io as sio
import os
# os.environ['PROJ_LIB'] = 'C:/Users/antonioh/Anaconda3/Library/share/proj/'
from matplotlib.colors import ListedColormap#, LinearSegmentedColormap
import netCDF4 as nc

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
            
def inter2points2nc(folder_tmp,file_name,time_waves,points_all,point_type,triang,Hs,Tm,Tp,Dir,Windv_x,Windv_y):
    
    
    Hsp = np.zeros((len(points_all),len(time_waves)))
    Tmp = np.zeros((len(points_all),len(time_waves)))
    Tpp = np.zeros((len(points_all),len(time_waves)))
    Dirp = np.zeros((len(points_all),len(time_waves)))
    Windv_xp= np.zeros((len(points_all),len(time_waves)))
    Windv_yp = np.zeros((len(points_all),len(time_waves)))
    
    for ts in range(0,len(Hs)):
        
        #interpolate to the points surroounding the land but still too deep for wave breaking
        intsurf = tri.LinearTriInterpolator(triang,Hs[ts])
        Hsp[:,ts] =intsurf(points_all[:,0],points_all[:,1])
        intsurf = tri.LinearTriInterpolator(triang,Tp[ts])
        Tpp[:,ts] =intsurf(points_all[:,0],points_all[:,1])
        intsurf = tri.LinearTriInterpolator(triang,Tm[ts])
        Tmp[:,ts] =intsurf(points_all[:,0],points_all[:,1])
        intsurf = tri.LinearTriInterpolator(triang,Dir[ts])
        Dirp[:,ts] =intsurf(points_all[:,0],points_all[:,1])
        intsurf = tri.LinearTriInterpolator(triang,Windv_x[ts])
        Windv_xp[:,ts] =intsurf(points_all[:,0],points_all[:,1])
        intsurf = tri.LinearTriInterpolator(triang,Windv_y[ts])
        Windv_yp[:,ts] =intsurf(points_all[:,0],points_all[:,1])
    
    fn = folder_tmp  + file_name
    try:
        os.remove(fn)
        print('removing wind_and_waves.nc of previous run')
    except:
        print('saving time serie at buoy location')
        
    ds = nc.Dataset(fn, 'w', format='NETCDF4')
    time = ds.createDimension('time', None)
    index = ds.createDimension('index', None)
    times = ds.createVariable('time', 'f8', ('time',))
    times.units='hours since 1950-01-01 00:00:00'
    times.calendar='gregorian'
    times[:] = [nc.date2num(x,units=times.units,calendar=times.calendar) for x in time_waves]
    lonnc = ds.createVariable('lon', 'f4', ('index',))
    lonnc.units ='degrees_east'
    lonnc[:] = points_all[:,0]
    latnc = ds.createVariable('lat', 'f4', ('index',))
    latnc.units ='degrees_north'
    latnc[:] = points_all[:,1]
    Hsnc= ds.createVariable('Hs', 'f4', ('index','time'))
    Hsnc.units = 'm'
    Hsnc[:,:] =  Hsp
    Tpnc= ds.createVariable('Tp', 'f4', ('index','time'))
    Tpnc.units = 's'
    Tpnc[:,:] =  Tpp
    Tmnc= ds.createVariable('Tm', 'f4', ('index','time'))
    Tmnc.units = 's'
    Tmnc[:,:] =  Tmp
    Dirnc= ds.createVariable('Dir', 'f4', ('index','time'))
    Dirnc.units = 'degress from north (north=0, east=90)'
    Dirnc[:,:] =  Dirp
    Windxnc= ds.createVariable('Windx', 'f4', ('index','time'))
    Windxnc.units = 'm/s'
    Windxnc[:,:] =  Windv_xp
    Windync= ds.createVariable('Windy', 'f4', ('index','time'))
    Windync.units = 'm/s'
    Windync[:,:] =  Windv_yp
    ds.setncattr_string('point_type', point_type)
    ds.close()

def plot_Hs_Dir_maps(nat,tar,dtt,x,y,triang,Time,Hs,Dir,result_folder):

    f = open('../extras/Hs_colormap.dat', 'r') # 'r' = read
    colors = np.genfromtxt(f, delimiter='  ')
    f.close()
    colors= np.hstack((colors,np.ones((len(colors[:,1]),1))))
    newcmp = ListedColormap(colors)

    ## Mesh grid values for interpolation
    # mesh grid national scale

    # linear interpolation - Dir
    Dir_cor = (270-Dir)*np.pi/180
    xDir = np.cos(Dir_cor)
    yDir = np.sin(Dir_cor)
    
    
    levels = np.arange(0.,4, 0.1)
    
    
    # Fongafale grid
    xllt = 178.84
    yllt = -8.77
    xurt = 179.39
    yurt = -8.33
    sgrid_x, sgrid_y = np.mgrid[xllt:xurt:30j,yllt:yurt:30j]
    
    # Tuvalu grid
    xllg = 174.5
    yllg = -11.5
    xurg = 181.5
    yurg = -4.5
    lgrid_x, lgrid_y = np.mgrid[xllg:xurg:20j,yllg:yurg:20j]


    
    plt.ioff()
    for ts in range(24,len(Hs),dtt):  
        
        fxDir = tri.LinearTriInterpolator(triang,xDir[ts])
        fyDir = tri.LinearTriInterpolator(triang,yDir[ts])
        
        if nat==1:
            lxdir = fxDir(lgrid_x,lgrid_y)
            lydir = fyDir(lgrid_x,lgrid_y)

            fig, ax = plt.subplots(figsize=(10,10))
            ax.set_aspect('equal')
            tcf = ax.tricontourf(triang, Hs[ts],levels=levels,cmap=newcmp,zorder=1)
            cbar = fig.colorbar(tcf, ax=ax)
            cbar.solids.set_edgecolor("face")
            ax.quiver(lgrid_x,lgrid_y,lxdir,lydir,scale_units='xy',angles='xy',zorder=2) #Dir x and 
            ax.set_title('Kiribati')
            ax.set_xlabel('Longitude (degrees)')
            ax.set_ylabel('Latitude (degrees)')
            ax.set_xlim(xllg,xurg)
            ax.set_ylim(yllg,yurg)
            plt.title(Time[ts])
            fileprint=result_folder + 'Tuvalu_' + Time[ts]
            plt.savefig(fileprint)
            #plt.show()
            plt.close(fig)
    
        if tar==1:
            

            sxdir = fxDir(sgrid_x,sgrid_y)
            sydir = fyDir(sgrid_x,sgrid_y)

            
            fig2, ax2 = plt.subplots(figsize=(10,10))
            ax2.set_aspect('equal')
            tcf2 = ax2.tricontourf(triang, Hs[ts],levels=levels,cmap=newcmp,zorder=1)
            cbar2 = fig2.colorbar(tcf2, ax=ax2)
            cbar2.solids.set_edgecolor("face")
            ax2.quiver(sgrid_x,sgrid_y,sxdir,sydir,scale_units='xy',angles='xy',zorder=2) #Dir x and y
            ax2.set_xlabel('Longitude (degrees)')
            ax2.set_ylabel('Latitude (degrees)')
            ax2.set_xlim(xllt,xurt)
            ax2.set_ylim(yllt,yurt)
            plt.title(Time[ts])
            fileprint=result_folder + 'Fongagale_' + Time[ts]
            plt.savefig(fileprint)
            plt.close(fig2)


def SWAN2nc(x,y,triang,time_waves,Hs,Tm,Tp,Dir,Windv_x,Windv_y,xll,yll,xur,yur,inc,namenc):

    
    x = np.arange(xll,xur,inc) 
    y = np.arange(yll,yur,inc)
    lgrid_x, lgrid_y = np.meshgrid(x,y)
    

    Hsp = np.ones((len(time_waves),lgrid_x.shape[0],lgrid_x.shape[1]))*-999
    Tmp = np.ones((len(time_waves),lgrid_x.shape[0],lgrid_x.shape[1]))*-999
    Tpp = np.ones((len(time_waves),lgrid_x.shape[0],lgrid_x.shape[1]))*-999
    Dirp = np.ones((len(time_waves),lgrid_x.shape[0],lgrid_x.shape[1]))*-999 
    Windp = np.ones((len(time_waves),lgrid_x.shape[0],lgrid_x.shape[1]))*-999
    WindDirp = np.ones((len(time_waves),lgrid_x.shape[0],lgrid_x.shape[1]))*-999

    for ts in range(48,len(time_waves)):# change the starting hour into range to the lenght of the spin-up 
        
        intsurf = tri.LinearTriInterpolator(triang,Hs[ts])
        Hsp[ts,:,:] =intsurf(lgrid_x,lgrid_y)
        intsurf = tri.LinearTriInterpolator(triang,Tp[ts])
        Tpp[ts,:,:] =intsurf(lgrid_x,lgrid_y)
        intsurf = tri.LinearTriInterpolator(triang,Tm[ts])
        Tmp[ts,:,:] =intsurf(lgrid_x,lgrid_y)
        intsurf = tri.LinearTriInterpolator(triang,Dir[ts])
        Dirp[ts,:,:] =intsurf(lgrid_x,lgrid_y)
        intsurf = tri.LinearTriInterpolator(triang,Windv_x[ts])
        wx=intsurf(lgrid_x,lgrid_y)
        intsurf = tri.LinearTriInterpolator(triang,Windv_y[ts])
        wy =intsurf(lgrid_x,lgrid_y)
        Windp[ts,:,:]= np.sqrt(wx**2 + wy**2)
        WindDirp[ts,:,:] = 270 - np.arctan(wy/wx)*180/np.pi
  
    
    ds = nc.Dataset(namenc, 'w', format='NETCDF4')
    time = ds.createDimension('time', None)
    lon = ds.createDimension('lon', None)
    lat = ds.createDimension('lat', None)
    times = ds.createVariable('time', 'f4', ('time',))
    times.units='hours since 1950-01-01 00:00:00'
    times.calendar='gregorian'
    times[:] = [nc.date2num(ti,units=times.units,calendar=times.calendar) for ti in time_waves]
    lonnc =ds.createVariable('lon', 'f4', ('lon',))
    lonnc.units ='degrees_east'
    lonnc[:] = x
    latnc = ds.createVariable('lat', 'f4', ('lat',))
    latnc.units ='degrees_north'
    latnc[:] = y
    Hsnc= ds.createVariable('Hs', 'f4', ('time','lat','lon'))
    Hsnc.units = 'm'
    Hsnc[:] =  Hsp
    Tpnc= ds.createVariable('Tp', 'f4', ('time','lat','lon'))
    Tpnc.units = 's'
    Tpnc[:] =  Tpp
    Tmnc= ds.createVariable('Tm', 'f4', ('time','lat','lon'))
    Tmnc.units = 's'
    Tmnc[:] =  Tmp
    Dirnc= ds.createVariable('Dir', 'f4', ('time','lat','lon'))
    Dirnc.units = 'degress from north (north=0, east=90)'
    Dirnc[:] =  Dirp
    Windc= ds.createVariable('Wind', 'f4', ('time','lat','lon'))
    Windc.units = 'm/s'
    Windc[:] =  Windp
    DirWindc= ds.createVariable('DirWind', 'f4', ('time','lat','lon'))
    DirWindc.units = 'degress from north (north=0, east=90)'
    DirWindc[:] =  WindDirp
    ds.close()
    
    

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

def read_loc_csv(fpath):
    df = pd.read_csv(fpath)
    atoll = df['atoll'].tolist()
    xll = df['xll'].tolist()
    yll = df['yll'].tolist()
    xur = df['xur'].tolist()
    yur = df['yur'].tolist()
    return(atoll,xll,yll,xur,yur)

def buoyhindcastforecast(now):
    folder_name='../archives/'
    flist = os.listdir(folder_name)
    dates = []
    for file in flist:
        if file.endswith("_buoy.nc"):
            nb = file.split("_")[0]
            dates.append(dt.datetime.strptime(nb+'0000',"%Y%m%d%H%M%S"))       
    dates.sort()
    

    output =open('../tmp/buoy_hindcast_forecast.txt','w')
    for r in dates:
         
        nc_fname = folder_name +  r.strftime("%Y") + r.strftime("%m") + r.strftime("%d") + r.strftime("%H") + '_buoy.nc'
        time = readvars_nc(nc_fname,'time')
        year = np.array([t.strftime('%Y') for t in time])
        month = np.array([t.strftime('%m') for t in time])
        day = np.array([t.strftime('%d') for t in time])
        hour = np.array([t.strftime('%H') for t in time])
        Hs = readvars_nc(nc_fname,'Hs')
        Tm = readvars_nc(nc_fname,'Tm')
        Tp = readvars_nc(nc_fname,'Tp')
        Dir = readvars_nc(nc_fname,'Dir')
        Wx = readvars_nc(nc_fname,'Windx')
        Wy = readvars_nc(nc_fname,'Windy')
        Wind = (Wx**2+Wy**2)**(1/2)
        WindDir= 270 - (np.arctan2(Wy,Wx)*180/np.pi)
        
        if r != dates[-1]:
            for i in range(47,53):
                s="{}".format(year[i])+" {}".format(month[i])+" {}".format(day[i])+" {}".format(hour[i])+ \
                " {:.3f}".format(Hs[0,i]) +" {:.3f}".format(Tm[0,i]) + " {:.3f}".format(Tp[0,i])+" {:.3f}".format(Dir[0,i])  +\
                " {:.3f}".format(Wind[0,i]) +" {:.3f}".format(WindDir[0,i]) + '\n'
                output.write(s)
                          
        elif r==dates[-1]:
            
            for i in range(47,len(time)):
                s="{}".format(year[i])+" {}".format(month[i])+" {}".format(day[i])+" {}".format(hour[i])+ \
                " {:.3f}".format(Hs[0,i]) +" {:.3f}".format(Tm[0,i]) + " {:.3f}".format(Tp[0,i])+" {:.3f}".format(Dir[0,i])  +\
                " {:.3f}".format(Wind[0,i]) +" {:.3f}".format(WindDir[0,i]) + '\n'
                output.write(s)
    output.close()
        

###################################################################################################################
def postprocess_SWAN(now,plot):

    folder_tmp ='../tmp/'
    out_name='../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H")  +'/'
    csvpath = '../extras/Tuvalu_island_gauge_locations.csv'
    result_folder = out_name + '/results/'
    try :
        os.mkdir(result_folder)
    except:
        print("results already exists")
        
        
    fl_name = out_name +'output.mat'

    # f = open('../extras/Lagoon_profiles_xyxy.txt', 'r') 
    # points_lago = np.genfromtxt(f, delimiter='  ')
    # f.close()
    # points_lago = points_lago[:,2:4]
    # # npointslago = len(points_lago)
    # f = open('../extras/Forereef_profiles_xyxy_v2.txt', 'r') 
    # points_reef = np.genfromtxt(f, delimiter='  ')
    # f.close()
    # # npointsreef = len(points_reef)
    # points_reef = points_reef[:,2:4]
    
    # point12 = [[172.79456],[1.4649415]]
    # pointE = [[173.06],[1.48]]
    # pointS = [[173],[1.32]]
    # pointW = [[172.88],[1.45]]
    
    
    # points_all=np.concatenate((points_lago,points_reef,np.array(point12).T,np.array(pointE).T,np.array(pointS).T,np.array(pointW).T),axis=0)
    # point_type= ['lagoon' for i in range(0,len(points_lago))] + ['forereef' for i in range(0,len(points_reef))] + ['point_12'] + ['point_E'] + ['point_S'] + ['point_W']
    
    f14 = out_name + 'fort.14'
    # open mesh file
    mesh = AdcircMesh.open(f14,crs=None)
    
    x=np.array(mesh.x)
    y=np.array(mesh.y)
    
    
    Time = load_SWAN_data(fl_name,"time")
    time_waves = [dt.datetime.strptime(Time[i],'%Y%m%d_%H%M%S') for i in range(0,len(Time))]
    
    
    Hs = load_SWAN_data(fl_name,"Hsig")
    Tp = load_SWAN_data(fl_name,"RTpeak")
    Tm = load_SWAN_data(fl_name,"Tm02")
    Dir = load_SWAN_data(fl_name,"Dir_")
    Windv_x = load_SWAN_data(fl_name,"Windv_x")
    Windv_y = load_SWAN_data(fl_name,"Windv_y")
    
    
    Hs = np.nan_to_num(Hs)
    Tp = np.nan_to_num(Tp)
    Tm = np.nan_to_num(Tm)
    Dir = np.nan_to_num(Dir)
    Windv_x = np.nan_to_num(Windv_x)
    Windv_y = np.nan_to_num(Windv_y)
   

    
     ## Trimesh triangulation used for plots and interpolation
    triang = tri.Triangulation(x,y,mesh.triangles)
    
  
    ## Mesh grid values for interpolation national scale Tuvalu grid boundary limit coordinates

    xll = 174.5
    yll = -11.5
    xur = 181.5
    yur = -4.5
    inc =0.01# about1km
    namenc = out_name  + '/results/Tuvalu.nc'  
    SWAN2nc(x,y,triang,time_waves,Hs,Tm,Tp,Dir,Windv_x,Windv_y,xll,yll,xur,yur,inc,namenc)
    
    
    ## Mesh grid values for interpolation island scale grid boundary limit coordinates
    atoll_list,xll_list,yll_list,xur_list,yur_list = read_loc_csv(csvpath)
    for i in range(len(atoll_list)):
        print("\n",atoll_list[i])
        inc =0.0005# about100m
        namenc = out_name  + '/results/P'+str(i+1)+'_'+str(atoll_list[i])+'.nc'
        SWAN2nc(x,y,triang,time_waves,Hs,Tm,Tp,Dir,Windv_x,Windv_y,xll_list[i],yll_list[i],xur_list[i],yur_list[i],inc,namenc)
        
  
    
    # file_name='wind_and_waves.nc'
    # inter2points2nc(folder_tmp,file_name,time_waves,points_all,point_type,triang,Hs,Tm,Tp,Dir,Windv_x,Windv_y)
    
    buoy =np.array([[179.0696333,   -8.4656333 ]])
    archfolder='../archives/' 
    file_name = now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") + now.strftime("%H") + '_buoy.nc'
    inter2points2nc(archfolder,file_name,time_waves,buoy,'buoy',triang,Hs,Tm,Tp,Dir,Windv_x,Windv_y)
    
    
    buoyhindcastforecast(now)

    if plot==1:
        
        nat=1 #plot national scale
        tar=0 #plot Tarawa scale
        dtt = 1 #time interval in hours
      
        plot_Hs_Dir_maps(nat,tar,dtt,x,y,triang,Time,Hs,Dir,result_folder)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    