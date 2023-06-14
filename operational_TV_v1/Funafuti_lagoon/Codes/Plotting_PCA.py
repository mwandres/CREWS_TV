# basics
import os
import os.path as op
import copy
from datetime import datetime, date
import calendar

# pip
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import matplotlib.colors as colors
from cftime._cftime import DatetimeGregorian
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy

# common fig parameters
_faspect = 1.618
_fsize = 9.8
_fdpi = 128

# TODO: figure parameters
_fntsize_label = 8
_fntsize_legend = 8
_fntsize_title = 8


def plot_agg_fams(dataset):
    '''
    This function plots the families that were previously created with the 
    basic_functions.py functions

    Parameters
    ----------
    dataset : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    fig, axes = plt.subplots(ncols=3,nrows=3,figsize=(12,12), 
                             subplot_kw=dict(polar=True))
    norm = colors.LogNorm(0.001,10)
    cmap = cm.jet
    dirs, tps = [], []
    for var in dataset:
        if 'Dir' in var:
            dirs.append(var)
        elif 'Tp' in var:
            tps.append(var)
    # plot all the partitions
    for axi, ax in enumerate(axes.flatten()):
        if axi==7: 
            break
        not_na = ~np.isnan(dataset[dirs[axi]].values)
        ax.hist2d(np.deg2rad(dataset[dirs[axi]].values[not_na]),
                  dataset[tps[axi]].values[not_na],
                  density=True,bins=[50,50],norm=norm,cmap=cmap)
        ax.set_title('Partition {}'.format(axi),fontsize=16,y=1.1)
        ax.set_theta_zero_location('N',offset=0)
        ax.set_theta_direction(-1)
        ax.set_facecolor('darkblue')
    cb = fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),
                      cax=fig.add_axes([0.95,0.15,0.02,0.7]),
                      orientation='vertical',alpha=0.4)
    cb.set_label('Density',rotation=270,fontsize=16)
    fig.suptitle('Normalized HISTOGRAMS of [Dir, Tp] \n sea=0, swells=rest', 
                 fontsize=20,y=1.02)
    fig.patch.set_facecolor('white')


def Plot_EOFs_EstelaPred(xds_PCA, n_plot, monthly=None, mask_land=None, show=True):
    '''
    Plot annual EOFs for 3D predictors

    xds_PCA:
        (n_components, n_components) PCs
        (n_components, n_features) EOFs
        (n_components, ) variance

        (n_lon, ) pred_lon: predictor longitude values
        (n_lat, ) pred_lat: predictor latitude values
        (n_time, ) pred_time: predictor time values

        method: gradient + estela

    n_plot: number of EOFs plotted
    '''

    # TODO: fix data_pos, fails only after pred.Load()?

    # PCA data
    variance = xds_PCA['variance'].values[:]
    EOFs = np.transpose(xds_PCA['EOFs'].values[:])
    PCs = np.transpose(xds_PCA['PCs'].values[:])
    data_pos = xds_PCA['pred_data_pos'].values[:]  # for handling nans
    pca_time = xds_PCA['pred_time'].values[:]
    pred_name = xds_PCA.attrs['pred_name']

    # PCA lat lon metadata
    lon = xds_PCA['pred_lon'].values
    lat = xds_PCA['pred_lat'].values

    # percentage of variance each field explains
    n_percent = variance / np.sum(variance)

    l_figs = []
    for it in range(n_plot):

        # get vargrd 
        var_grd_1d = EOFs[:,it] * np.sqrt(variance[it])

        # insert nans in data
        base = np.nan * np.ones(data_pos.shape)
        base[data_pos] = var_grd_1d

        var = base[:int(len(base)/2)]
        grd = base[int(len(base)/2):]

        # reshape data to grid
        C1 = np.reshape(var, (len(lon), len(lat)))
        C2 = np.reshape(grd, (len(lon), len(lat)))

        # figure
        fig = plt.figure(figsize=(_faspect*_fsize, 2.0/3.0*_fsize))

        # layout
        gs = gridspec.GridSpec(4, 4, wspace=0.10, hspace=0.2)

        ax_EOF_1 = plt.subplot(gs[:3, :2])
        ax_EOF_2 = plt.subplot(gs[:3, 2:])
        ax_evol = plt.subplot(gs[3, :])

        # EOF pcolormesh (SLP and GRADIENT)
        axplot_EOF(ax_EOF_1, C1, lon, lat, ttl = 'SLP', land=mask_land)
        axplot_EOF(ax_EOF_2, C2, lon, lat, ttl = 'GRADIENT', land=mask_land)

        # time series EOF evolution
        evol =  PCs[it,:]/np.sqrt(variance[it])
        if monthly:
            monthly_evol = monthly[:,it]/np.sqrt(variance[it])
            axplot_EOF_evolution(ax_evol, pca_time, evol, monthly_evol=monthly_evol)
        else:
            axplot_EOF_evolution(ax_evol, pca_time, evol)

        # figure title
        ttl = 'EOF #{0}  ---  {1:.2f}%'.format(it+1, n_percent[it]*100)
        fig.suptitle(ttl, fontsize=14, fontweight='bold')

        l_figs.append(fig)

    # show and return figure
    if show: plt.show()
    return l_figs


def axplot_EOF(ax, EOF_value, lon, lat, ttl='', land=None):
    'axes plot EOFs 2d map'

    cmap = cm.get_cmap('RdBu_r')

    # EOF pcolormesh 
    ax.pcolormesh(
        lon, lat, np.transpose(EOF_value),
        cmap=cmap, shading='gouraud',
        clim=(-1,1),
    )

    # optional mask land
    if type(land).__module__ == np.__name__:
        add_land_mask(ax, lon, lat, land, 'grey')

    # axis and title
    ax.set_title(
        ttl,
        {'fontsize': 10, 'fontweight':'bold'}
    )
    ax.tick_params(axis='both', which='major', labelsize=8)

    
def axplot_EOF_evolution(ax, time, EOF_evol, monthly_evol=None):
    'axes plot EOFs evolution'

    # date axis locator
    yloc1 = mdates.YearLocator(1)
    yfmt = mdates.DateFormatter('%Y')

    # convert to datetime
    dtime = [n2d(t) for t in time]

    # plot EOF evolution 
    ax.plot(
        dtime, EOF_evol,
        linestyle='-', linewidth=0.5, color='black',
        label='Daily PCs'
    )
    
    if monthly_evol:
        monthly_evol.plot(ax=ax,color='red',label='Monthly mean',lw=2)

    # configure axis
    ax.set_title('')
    ax.set_xlim(time[0], time[-1])
    ax.xaxis.set_major_locator(yloc1)
    ax.xaxis.set_major_formatter(yfmt)
    ax.grid(True, which='both', axis='x', linestyle='--', color='grey')
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.legend(ncol=2)

    
def add_land_mask(ax, lon, lat, land, color):
    'addsland mask pcolormesh to existing pcolormesh'

    # select land in mask
    landc = land.copy()
    landc[np.isnan(land)]=1
    landc[land==1]=np.nan

    ax.pcolormesh(
        lon, lat, landc,
        cmap=colors.ListedColormap([color]), shading='gouraud',
    )

    
def n2d(dt64):
    'converts np.datetime64[ns] into datetime'

    ts = (dt64 - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
    return datetime.utcfromtimestamp(ts)



