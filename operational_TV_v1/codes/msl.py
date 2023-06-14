import scipy.io as sio
import numpy as np
from matplotlib import pyplot as plt
import datetime
from scipy.interpolate import interp2d
import netCDF4
from netCDF4 import num2date

   
def cm2inch(*tupl):
    inch = 2.54
    if type(tupl[0]) == tuple:
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple()
    
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
        #print(i)
        nc_fnameX = nc_fname + '_' +  "{0:0>3}".format(i) +'.nc'
        nc = netCDF4.Dataset(nc_fnameX)
        
        tt = np.append(tt,num2date(np.array(nc['time']), units=nc['time'].units, calendar = nc['time'].calendar),axis=0)
        Uwind = np.append(Uwind,np.array(nc['u-component_of_wind_height_above_ground']),axis=0)
        Vwind = np.append(Vwind,np.array(nc['v-component_of_wind_height_above_ground']),axis=0)
        MSLP =  np.append(MSLP,np.array(nc['Pressure_reduced_to_MSL_msl']),axis=0)

        #print(nc_fnameX)
    Uwind = np.squeeze(Uwind)
    Vwind = np.squeeze(Vwind)
    MSLP = np.squeeze(MSLP)
    return(ln,lt,tt,Uwind,Vwind,MSLP)


def plot_winds_and_MSLP_report(i, wnd_nc):
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
    
    fU = interp2d(ln, lt, Uwind[i,:,:], kind='cubic')
    UU = fU(grid_x, grid_y)
    fV = interp2d(ln, lt, Vwind[i,:,:], kind='cubic')
    VV = fV(grid_x, grid_y)
    
    levels = np.arange(0.,40, 0.1)
    levels_MSLP = np.arange(970.,1030, 10)
    
    return (tt,ln, lt, WndSpd, levels,lgrid_x,grid_x,lgrid_y,grid_y,Uwind,Vwind,MSLP,levels_MSLP)


def plotter_msl(ax2, time_step, bottomlabel,labelleftt, havebar,fig,wnd_nc):
    ##GET DATA
    tt,ln,lt,windspeed,levels,lgrid_x,grid_x,lgrid_y,grid_y,Uwind,Vwind,MSLP,levels_MSLP = plot_winds_and_MSLP_report(time_step,wnd_nc)
    fU = interp2d(ln, lt, Uwind[time_step,:,:], kind='cubic')
    UU = fU(grid_x, grid_y)
    fV = interp2d(ln, lt, Vwind[time_step,:,:], kind='cubic')
    VV = fV(grid_x, grid_y)

    tcf = ax2.contourf(ln,lt,windspeed[time_step,:,:],levels=levels,cmap='gist_ncar')
    if havebar:
        cbar = fig.colorbar(tcf, ax=ax2)
        cbar.ax.tick_params(labelsize=6)
        cbar.set_label('Wind Speed (m/s)',fontsize=6)
    
    ax2.quiver(lgrid_x,lgrid_y,UU,VV)
    CS = ax2.contour(ln,lt,MSLP[time_step,:,:]/100,levels=levels_MSLP,colors='w')
    ax2.clabel(CS, inline=True, fontsize=10)
    if bottomlabel:
        ax2.set_xlabel('Longitude (degrees)', fontsize=7)
    
    ax2.set_title(str(tt[time_step]),fontsize=7)
    ax2.tick_params(axis='both', which='both', labelsize=6)
    #ax2.tick_params(top=False, which='both',axis='both',labeltop=False, bottom=False, labelbottom=False)
    if labelleftt:
        ax2.set_ylabel('Latitude (degrees)',fontsize=7)
    ax2.tick_params(left=labelleftt, right=False, labelleft = labelleftt, labelbottom = bottomlabel, bottom = bottomlabel)
    return ax2

figsize = cm2inch((21,29.7))

def msl_plots(factor, wnd_nc):
    fig, ax = plt.subplots(figsize=figsize, dpi=100)
    ax.axis('off')
    ax.text(0.01, 1.09,"Regional MSLP Maps", transform=ax.transAxes,fontsize=8, verticalalignment='top', color='#2F5F8A')
    #SECOND SET
    ax5 = fig.add_axes([0.06, 0.845, 0.30, 0.20])
    ax5 = plotter_msl(ax5, 0, False, True,False,fig,wnd_nc)
    ax6 = fig.add_axes([0.37, 0.845, 0.30, 0.20])
    ax6 = plotter_msl(ax6, 2, False, False,False,fig,wnd_nc)
    ax7 = fig.add_axes([0.68, 0.845, 0.30, 0.20])
    ax7 = plotter_msl(ax7, 4, False, False,True,fig,wnd_nc)

    #THIRD SET
    ax8 = fig.add_axes([0.06, 0.625, 0.30, 0.20])
    ax8 = plotter_msl(ax8, 6, False, True,False,fig,wnd_nc)
    ax9 = fig.add_axes([0.37, 0.625, 0.30, 0.20])
    ax9 = plotter_msl(ax9, 8, False, False,False,fig,wnd_nc)
    ax10 = fig.add_axes([0.68, 0.625, 0.30, 0.20])
    ax10 = plotter_msl(ax10, 10, False, False,True,fig,wnd_nc)

    #LAST SET
    ax11 = fig.add_axes([0.06, 0.405, 0.375, 0.20])
    ax11 = plotter_msl(ax11, 12, True, True,True,fig,wnd_nc)
    fig.savefig('msl.png', bbox_inches='tight')
    plt.close(fig)
    return None

#fig = msl_plots()
#fig.savefig('msl.png', bbox_inches='tight')
