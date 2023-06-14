import os
import os.path as op
import sys

import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize,ListedColormap, LinearSegmentedColormap
import matplotlib as mpl
from matplotlib import cm

import xarray as xr
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import pickle as pk
from scipy.stats import entropy

import cartopy
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

sys.path.insert(0, op.join(os.path.abspath(''), '..'))
from Codes.pca import *
from Codes.Plotting_PCA import *


import numpy as np



#########################################################################################################
##                                             CALCULATIONS                                            ##
#########################################################################################################

def Set_t_h_sp(sp):
    
    sp['t']=(['freq'],1/sp.freq)    
    h=np.full(np.shape(sp.efth.values),0.0)
    freq_d=np.diff(sp.freq)
    freq_d=np.append(0.00001,freq_d)
    for i in range(len(sp.freq)):
        h[:,i,:]=((sp.efth.values[:,i,:])*(sp.dir.values[0]-sp.dir.values[1])*freq_d[i])
    sp['h']=(['time','freq','dir'],h)
    
    return sp

def sp_group_bins(sp, dir_bins, t_bins):
    
    print('Direction bins:' + str(dir_bins))
    print('---------------------------------')
    print('Period bins:' + str(t_bins))
    
    sp_mod = xr.Dataset(
        {
            'h': (['time','t','dir'], sp.h.values),
        },
        coords={
            't': sp.t.values,
            'dir': sp.dir.values,
            'time': sp.time.values,
        },
    )
    sp_mod=sp_mod.groupby_bins('dir',dir_bins).sum()
    sp_mod=sp_mod.groupby_bins('t',t_bins).sum()
    sp_mod['dir']=('dir_bins',(dir_bins)[1:])#+(0.5*(dir_bins[0]-dir_bins[1])))[1:])
    sp_mod['t']=('t_bins',t_bins[1:])
    # sp_mod['h']=4*np.sqrt(sp_mod['h'])
    sp_mod['h_t']=4*np.sqrt(sp_mod['h'])
    
    return sp_mod



def spec_kmeans_clustering(path_main,site, sp_mod_daily, n_pcs, num_clusters, do=[], min_data_frac=10):
    
    if do==1:
        min_data=np.int(len(sp_mod_daily.time)/num_clusters/min_data_frac)
        num, it=0,0

        while np.nanmin(num)<min_data:

            sys.stdout.write('\rIteration: %d' %it)
            sys.stdout.flush()

            kma = KMeans(n_clusters=num_clusters,init='k-means++',n_init=10,tol=0.0000000001) #80
            bmus = kma.fit_predict(sp_mod_daily.PCs.values[:,:n_pcs])
            centers= kma.cluster_centers_
            num=[]
            for i in range(num_clusters):
                num=np.append(num,len(np.where(bmus==i)[0]))
            it=it+1

        # km, x and var_centers
        kma_order = np.argsort(np.mean(-centers, axis=1))

        sorted_bmus = np.zeros((len(bmus),),)*np.nan
        for i in range(num_clusters):
            posc = np.where(bmus == kma_order[i])
            sorted_bmus[posc] = i

        bmus=sorted_bmus

        pk.dump(kma, open(os.path.join(path_main, 'Results', site , 'Spec_KMA_' + site + '.pkl'),"wb"))

        print('Number of sims: ' + str(kma.n_iter_))
        _, group_size = np.unique(bmus, return_counts=True)
        print('Minimum number of data: ' + str(np.nanmin(group_size)))
        sp_mod_daily['kma_order']=kma_order
        sp_mod_daily['n_pcs']=n_pcs
        sp_mod_daily['bmus']=('time',bmus)
        sp_mod_daily['centroids']=(('clusters','n_pcs_'),centers)
        sp_mod_daily.to_netcdf(os.path.join(path_main,'Results',site,'Spec_classification_PCS_KMA_SF.nc'))
    else:
        sp_mod_daily=xr.open_dataset(os.path.join(path_main,'Results',site,'Spec_classification_PCS_KMA_SF.nc'))

    return sp_mod_daily

def slp_kmeans_clustering(path_main, site, slp, PCs, n_pcs, estela_mask, num_clusters=[], do=1, min_data=20):
    
    if do:

        num, it=0,0

        while np.nanmin(num)<min_data:

            sys.stdout.write('\rIteration: %d' %it)
            sys.stdout.flush()

            kma = KMeans(n_clusters=num_clusters,init='k-means++',n_init=10,tol=0.0001) #80
            bmus = kma.fit_predict(PCs[:,:n_pcs])
            centers= kma.cluster_centers_
            num=[]
            for i in range(num_clusters):
                num=np.append(num,len(np.where(bmus==i)[0]))
            it=it+1

        pk.dump(kma, open(os.path.join(path_main, 'Results', site , 'SLP_KMA_' + site + '.pkl'),"wb"))

        print('Number of sims: ' + str(kma.n_iter_))
        _, group_size = np.unique(bmus, return_counts=True)
        print('Minimum number of data: ' + str(np.nanmin(group_size)))

        # km, x and var_centers
        kma_order = np.argsort(np.mean(-centers, axis=1))

        sorted_bmus = np.zeros((len(bmus),),)*np.nan
        for i in range(num_clusters):
            posc = np.where(bmus == kma_order[i])
            sorted_bmus[posc] = i

        bmus=sorted_bmus

        slp['bmus']=('time',bmus)
        slp['kma_order']=kma_order
        slp['cluster']=np.arange(num_clusters)
        slp['centers']=(('cluster','dims'),centers)
        slp['est_mask']=(('latitude','longitude'),estela_mask.values)
        slp
        slp.to_netcdf(os.path.join(path_main,'Results',site, 'SLP_classification_PCS_KMA_SF.nc'))
    else:
        slp=xr.open_dataset(os.path.join(path_main,'Results',site, 'SLP_classification_PCS_KMA_SF.nc'))
    return slp


def n34_format(n34, rolling_mean=[]):

    n34 = xr.Dataset(
        {
            'index': (['time'], n34[:,1:].reshape(-1)),
        },
        coords={
            'time': pd.date_range(datetime(n34[0,0].astype('int'),1,1), datetime(n34[-1,0].astype('int')+1,1,1) , freq='1M') ,
        },
    )
    n34.index[np.where(n34.index==-99.99)]=np.nan

    #Means are calculated from the period 1991-2020, same as the forecast
    n34['index']=n34.index-np.nanmean(n34.sel(time=slice('1991-01-01','2020-12-31')).index)
    
    n34_rol=n34.rolling(time=3, center=True).mean()
    if rolling_mean:
        ## 3-month rolling mean

        n34_rol=n34_rol.sel(time=slice('1979-01-01','2021-01-01'))
        n34=n34.sel(time=slice('1979-01-01','2021-01-01'))
        n34=n34_rol
    
    return n34


def calculate_conditioned_probabilities(slp, num_clusters_slp, sp, num_clusters_sp):
    
    C_T=np.full([np.sqrt(num_clusters_sp).astype('int'), np.sqrt(num_clusters_sp).astype('int'), num_clusters_slp],np.nan)

    for ic in range(num_clusters_slp):

        sel_2 = sp.bmus.values[np.where(slp.bmus==ic)[0][:]]

        # get DWT cluster probabilities
        cps = ClusterProbabilities(sel_2, range(num_clusters_sp))
        C_T[:,:,ic] = np.reshape(cps, (np.sqrt(num_clusters_sp).astype('int'), np.sqrt(num_clusters_sp).astype('int')))
        
    groups=np.reshape(range(num_clusters_sp), (np.sqrt(num_clusters_sp).astype('int'), np.sqrt(num_clusters_sp).astype('int')))
    conditioned_prob = xr.Dataset(
        {
            'prob_sp': (['ir','ic','bmus_slp'], C_T),
            'kma_sp': (['ir','ic'], groups),
        },
        coords={
            'ir': range(np.sqrt(num_clusters_sp).astype('int')),
            'ic': range(np.sqrt(num_clusters_sp).astype('int')),
            'bmus_slp': range(num_clusters_slp),
        },
    )
    return conditioned_prob


def calculate_entropy(sp_mod_daily, slp_dy_m, spec_kma):
    
    years=np.arange(2016,2021,1)
    months=np.arange(1,13,1)

    re_all=[]

    for y in years:
        for m in months:
            month=[str(y)+ '-' + str(m)]

            if (y==2020) & (m>10):
                continue
            else:
                sp=sp_mod_daily.sel(time=month[0])
                real_prob=np.full([len(spec_kma.clusters)],0.0)
                for a in range(len(spec_kma.clusters)):
                    s=np.where(sp.bmus.values==a)[0]
                    real_prob[a]=len(s)/len(sp.time)

                slp_m=slp_dy_m.sel(time=month[0]).isel(time=0)
                model_prob=slp_m.prob_spec.values.reshape(-1)

                real_prob[np.where(real_prob==0)]=0.00000000000000001
                model_prob[np.where(model_prob==0)]=0.00000000000000001
                RE=entropy(real_prob, model_prob)
                re_all.append(RE)
    return re_all


def download_generate_dynamic_forecast_data(path_main, slp_h, site, year, month, day, runs, est, do=1):
    if do:
        for run in runs:

            slp_f = download_prmsl(year, month, day,run, load_local_cache=False)
            slp_f=slp_f.drop(['time','meanSea','valid_time']).rename({'prmsl':'slp','ftime':'time'})
            slp_f=slp_f.sel(longitude=slp_h.longitude,latitude=slp_h.latitude).resample(time='1D').mean()
            slp_f=spatial_gradient(slp_f,'slp')

            est = est.sortby(est.latitude,ascending=False).interp(coords={'latitude':slp_f.latitude,
                                                                               'longitude':slp_f.longitude})
            estela_D = est.drop('time').traveltime
            estela_mask = est.estela_mask # mask for slpƒ

            pca_dy = dynamic_estela_predictor(slp_f, 'slp', estela_D)

            slp_f=slp_f.expand_dims(dim='run').assign_coords({'run':(('run'),[run])})
            pca_dy=pca_dy.expand_dims(dim='run').assign_coords({'run':(('run'),[run])})

            if run==0:
                SLP_F=slp_f
                PCA_DY=pca_dy
            else:
                SLP_F=xr.concat([SLP_F,slp_f],dim='run')
                PCA_DY=xr.concat([PCA_DY,pca_dy],dim='run')

        SLP_F.to_netcdf(os.path.join(path_main, 'Data', 'Forecast_prmsl_'+str(year) + '_' + str(month) +'_' + str(day) + '.nc'))
        PCA_DY.to_netcdf(os.path.join(path_main, 'Data', 'Dynamic_predictor_' + site + '_forecast_' +str(year) + '_' + str(month) +'_' + str(day) +'.nc'))
    else:
        SLP_F=xr.open_dataset(os.path.join(path_main, 'Data', 'Forecast_prmsl_'+str(year) + '_' + str(month) +'_' + str(day) + '.nc'))
        PCA_DY=xr.open_dataset(os.path.join(path_main, 'Data', 'Dynamic_predictor_' + site + '_forecast_' +str(year) + '_' + str(month) +'_' + str(day) +'.nc'))
    
    return SLP_F, PCA_DY

#########################################################################################################
##                                               PLOTTING                                              ##
#########################################################################################################

def Plot_Sp_Transformation(sp, sp_mod, time=[]):
    
    fig = plt.figure(figsize=[20,8])

    gs1=gridspec.GridSpec(1,1)
    ax1=fig.add_subplot(gs1[0],projection='polar')
    if time:
        z=sp.isel(time=time).efth.values
        lb=str(sp.isel(time=time).time.values)
    else:
        z=np.mean(sp.efth.values, axis=0)
        lb='Mean'
        
    p1=plot_spectrum(ax1,np.deg2rad(sp.dir.values), sp.freq.values,z)
    st=sp.isel(time=time).spec.stats(['hs','tp','dpm'])
    ax1.set_title(lb, fontsize=16)
    # ax1.set_title('sqrt (E): Bulk Hs=' + str(np.round(st.hs.values,2)) + ',   Dpm = ' + str(np.round(st.dpm.values,2)))
    gs1.tight_layout(fig, rect=[[], [], 0.4, []])

    gs2=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs2[0])
    plt.colorbar(p1,cax=ax0)
    ax0.set_ylabel('Sqrt(Efth)')
    gs2.tight_layout(fig, rect=[0.4, 0.1, 0.46, 0.9])

    gs3=gridspec.GridSpec(1,1)
    ax2=fig.add_subplot(gs3[0],projection='polar')
    if time:
        z=sp_mod.isel(time=time).h.values
        z=4*np.sqrt(z)
    else:
        z=np.mean(sp_mod.h.values,axis=2)
        z=4*np.sqrt(z)

    p2=plot_spectrum_hs(ax2,np.deg2rad(sp_mod.dir.values), sp_mod.t.values,z,vmax=1.75,ylim=0.99*np.max(sp_mod.t.values))
#     ax2.set_title('Hs (' + str(np.round(np.sqrt(np.sum(z**2)),2)) + 'm)' )
    ax2.set_title(lb, fontsize=16)
    gs3.tight_layout(fig, rect=[0.5, [], 0.9, []])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p2,cax=ax0)
    ax0.set_ylabel('Hs (m)')
    gs4.tight_layout(fig, rect=[0.92, 0.1, 0.98, 0.9])

def plot_pc_space(sp_mod_daily, figsize=[10,10]):
    
    fig = plt.figure(figsize=figsize)

    gs3=gridspec.GridSpec(2,2,hspace=0.1, wspace=0.1)
    ax=fig.add_subplot(gs3[0])
    ax.plot(sp_mod_daily.PCs.values[:,1],sp_mod_daily.PCs.values[:,0],'.', color='coral')
    ax.plot(sp_mod_daily.centroids[:,1],sp_mod_daily.centroids[:,0],'.', color='crimson',markersize=10)
    ax.set_ylabel('PC1', fontsize=14)
    ax.set_xlabel('PC2', fontsize=14)

    ax1=fig.add_subplot(gs3[1])
    ax1.plot(sp_mod_daily.PCs.values[:,2],sp_mod_daily.PCs.values[:,0],'.', color='coral')
    ax1.plot(sp_mod_daily.centroids[:,2],sp_mod_daily.centroids[:,0],'.', color='crimson',markersize=10)

    ax2=fig.add_subplot(gs3[3])
    ax2.plot(sp_mod_daily.PCs.values[:,2],sp_mod_daily.PCs.values[:,1],'.', color='coral', label='PCs')
    ax2.plot(sp_mod_daily.centroids[:,2],sp_mod_daily.centroids[:,1],'.', color='crimson',markersize=10, label='PCs - Centroids')
    ax2.set_ylabel('PC2', fontsize=14)
    ax2.set_xlabel('PC3', fontsize=14)
    ax2.legend(fontsize=12)

def plot_bmus_count(sp_mod_daily):
    
    plt.figure(figsize=[15,6])
    plt.grid(linestyle=':',alpha=0.4, color='darkmagenta')
    plt.hist(sp_mod_daily.bmus, np.arange(len(np.unique(sp_mod_daily.bmus))), color='indianred', alpha=0.4, edgecolor='yellow')
    plt.xlabel('Bmus', fontsize=15)
    plt.ylabel('Count', fontsize=15)

def Plot_Estela(est,extent, figsize=[20, 8]):
    
    fig = plt.figure(figsize=figsize)

    ax = plt.axes(projection = ccrs.PlateCarree(central_longitude=180))
    ax.set_extent(extent,crs = ccrs.PlateCarree())
    ax.stock_img()

    # cartopy land feature
    land_10m = cartopy.feature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='darkgrey', facecolor='gainsboro',  zorder=5)
    ax.add_feature(land_10m)
    ax.gridlines()

    # scatter data points
    vmin, vmax =96000, 103000

    p1=plt.pcolor(est.longitude.values, est.latitude.values, est.F.values, transform=ccrs.PlateCarree(),zorder=2, cmap='plasma')
    plt.contour(est.longitude.values, est.latitude.values, est.traveltime.values, levels=np.arange(1,26,1), transform=ccrs.PlateCarree(),zorder=2, linestyle=':', color='royalblue')

    plt.colorbar(p1).set_label('Energy')
    
def plot_spectrum(ax,x,y,z, vmax=0.3, ylim=0.49):
    x1=np.append(x,x[0])
    y1=np.append(0,y)
    z1=np.column_stack((z[:,:],z[:,-1]))
    p1=ax.pcolormesh(x1,y1,np.sqrt(z1), vmin=0, vmax=vmax)   
    p1.set_cmap('inferno')
    ax.set_theta_zero_location('N', offset=0)
    ax.set_theta_direction(-1)
    ax.set_ylim(0,ylim)
    ax.tick_params(axis='y', colors='plum',labelsize=14,grid_linestyle=':',grid_alpha=0.75,grid_color='plum')
    ax.tick_params(axis='x', colors='purple',labelsize=14,pad=5,grid_linestyle=':',grid_alpha=0.75,grid_color='plum')
    ax.grid(color='plum', linestyle='--', linewidth=0.7,alpha=0.2)
    return p1

def plot_random_daily_specs(sp_mod_daily, gr=6):
    
    fig = plt.figure(figsize=[17,13])

    gs3=gridspec.GridSpec(gr,gr,hspace=0.01, wspace=0.01)
    hs_mr=[]
    for b in range(gr**2):
        s=np.random.randint(len(sp_mod_daily.time), size=1)
        ax2=fig.add_subplot(gs3[b],projection='polar')
        z=sp_mod_daily.isel(time=s).h.values
        z=4*np.sqrt(z)
        z=np.reshape(z,(len(sp_mod_daily.t_bins),len(sp_mod_daily.dir_bins)))
        hs_mr.append(np.sqrt(np.sum(z**2)))
        p2=plot_spectrum_hs(ax2,np.deg2rad(sp_mod_daily.dir_bins.values), sp_mod_daily.t_bins.values,z,vmax=1.2,ylim=np.nanmax(sp_mod_daily.t_bins), remove_axis=1, cmap='magma') 
    gs3.tight_layout(fig, rect=[[], [], 0.9, []])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p2,cax=ax0)
    ax0.set_ylabel('Hs (m)', fontsize=14)
    gs4.tight_layout(fig, rect=[0.9, 0.1, 0.99, 0.9])
    
    
    plt.figure(figsize=[18,6])
    plt.grid(linestyle=':',alpha=0.4, color='darkmagenta')
    plt.plot(range(gr*gr), hs_mr,'.-', color='darkmagenta', alpha=0.8, markersize=15)
    plt.xlabel('Random case', fontsize=15)
    plt.ylabel('Hs_mean (m)', fontsize=15)

def plot_n34(n34, l1=0.5, l2=-0.5, figsize=[22,6]):
    plt.figure(figsize=figsize)
    plt.plot(n34.time,n34.index,'.-',color='grey', linewidth=2, markersize=8, alpha=0.5, label='Rolling mean')
    plt.grid(color='plum')
    plt.xlabel('Time',fontsize=16)
    plt.ylabel('Niño 3.4',fontsize=16)

    s1=np.where(n34.index>l1)[0]
    s2=np.where(n34.index<l2)[0]

    #Decision on wether to use the rolling mean or the original monthly values

    clasif=np.full([len(n34.index)],1.0)
    clasif[s1]=0
    clasif[s2]=2
    n34['classification']=('time',clasif) #0:Niño, 1:Neutral, 2:Niña

    plt.plot(n34.time[s1],n34.index[s1],'.',color='indianred', markersize=12,label='El Niño')
    plt.plot(n34.time[s2],n34.index[s2],'.',color='cornflowerblue', markersize=12,label='La Niña')


    plt.plot([n34.time.values[0],n34.time.values[-1]], [l1,l1],'-',color='indianred', linewidth=2, label='El Niño threshold')
    plt.plot([n34.time.values[0],n34.time.values[-1]], [0,0],'--',color='grey', linewidth=2)
    plt.plot([n34.time.values[0],n34.time.values[-1]], [l2,l2],'-',color='cornflowerblue', linewidth=2, label='La Niña threshold')
    plt.xlim([n34.time.values[0],n34.time.values[-1]])

    plt.legend(ncol=5, fontsize=13)

def plot_mjo(mjo, figsize=[7.5, 6]):
    
    fig, ax = plt.subplots(figsize=figsize)
    plt.plot([-4.5, 4.5],[0,0], color='grey')
    plt.plot([0,0],[-4.5, 4.5], color='grey')
    plt.plot([-4.5,4.5],[-4.5, 4.5], color='grey')
    plt.plot([-4.5, 4.5],[4.5,-4.5], color='grey')
    plt.xlim([-4.5,4.5])
    plt.ylim([-4.5,4.5])
    cax=plt.scatter(mjo.rmm1, mjo.rmm2,3, mjo.phase, cmap='Set2')
    plt.xlabel('RMM1', fontsize=16)
    plt.ylabel('RMM2', fontsize=16)
    tk=np.linspace(0.5,7.5,9)
    cbar=fig.colorbar(cax, ticks=tk)
    cbar.ax.set_yticklabels(['0','1','2','3','4','5','6','7','8']) 
    cbar.set_label('MJO phase', fontsize=15)
    
    
def plot_EOFS_spec(sp_mod_daily, gr=5):
    
    fig = plt.figure(figsize=[14,12])
    gs3=gridspec.GridSpec(gr,gr,hspace=0.01, wspace=0.01)
    for b in range(sp_mod_daily.n_pcs.values):
        ax2=fig.add_subplot(gs3[b],projection='polar')
        z=np.reshape(sp_mod_daily.EOFs.values[b,:],(len(sp_mod_daily.t_bins),len(sp_mod_daily.dir_bins)))
        p2=plot_spectrum_hs(ax2,np.deg2rad(sp_mod_daily.dir_bins.values), sp_mod_daily.t_bins.values,z,vmin=-0.6,vmax=0.6,ylim=np.nanmax(sp_mod_daily.t_bins), remove_axis=1, cmap='RdBu_r') 


    gs3.tight_layout(fig, rect=[[], [], 0.9, []])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p2,cax=ax0)
    ax0.set_ylabel('Hs (m)', fontsize=14)
    gs4.tight_layout(fig, rect=[0.9, 0.1, 0.99, 0.9])

def Plot_kmeans_clusters(sp_mod_daily, annomaly=[], prob_max=0.06, vmax=1.2):
    
    num_clusters=len(np.unique(sp_mod_daily.bmus))
    mean_z=4*np.sqrt(np.mean(sp_mod_daily.h.values,axis=0))
    
    fig = plt.figure(figsize=[16,12])

    h_clust=np.full([len(sp_mod_daily.t_bins), len(sp_mod_daily.dir_bins), len(sp_mod_daily.kma_order)], np.full)

    gs3=gridspec.GridSpec(np.int(np.sqrt(num_clusters)),np.int(np.sqrt(num_clusters)),hspace=0.01, wspace=0.01)
    hs_m=[]
    for b in range(num_clusters):
        s=np.where(sp_mod_daily.bmus==b)[0]
        prob=len(s)/len(sp_mod_daily.bmus)
        ax2=fig.add_subplot(gs3[b],projection='polar')
        if annomaly:
            z=4*np.sqrt(np.mean(sp_mod_daily.h.isel(time=s).values,axis=0)) - mean_z
            h_clust[:,:,b]=z
            p2=plot_spectrum_hs(ax2,np.deg2rad(sp_mod_daily.dir_bins.values), sp_mod_daily.t_bins.values,z,vmin=-0.3, vmax=vmax,ylim=np.nanmax(sp_mod_daily.t_bins), cmap='RdBu_r', remove_axis=1, prob=prob, prob_max=prob_max)  
        else:
            z=np.mean(sp_mod_daily.h.isel(time=s).values,axis=0)
            z=4*np.sqrt(z)
            p2=plot_spectrum_hs(ax2,np.deg2rad(sp_mod_daily.dir_bins.values), sp_mod_daily.t_bins.values, z, vmin=0, vmax=vmax, ylim=np.nanmax(sp_mod_daily.t_bins), remove_axis=1, cmap='CMRmap_r', prob=prob, prob_max=prob_max) 
            h_clust[:,:,b]=z
        hs_m.append(np.sqrt(np.sum(z**2)))
        
        
    gs3.tight_layout(fig, rect=[[], [], 0.83, []])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p2,cax=ax0)
    if annomaly:
        ax0.set_ylabel('Hs (m) - ANNOMALIES', fontsize=14)
    else:
        ax0.set_ylabel('Hs (m)', fontsize=14)
    gs4.tight_layout(fig, rect=[0.83, 0.1, 0.9, 0.9])

    gs5=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs5[0])
    norm = Normalize(vmin=0, vmax=prob_max)        
    cmap = cm.get_cmap('Blues')
    cb1 = mpl.colorbar.ColorbarBase(ax0, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_label('Probability', fontsize=14)
    gs5.tight_layout(fig, rect=[0.9, 0.1, 0.98, 0.9])
    
    return h_clust, hs_m
        
def Plot_kmeans_clusters_probability(sp_mod_daily, prob_max=0.06, vmax_z1=1.2):
    mean_z=np.mean(sp_mod_daily.h.values,axis=0)
    num_clusters=len(np.unique(sp_mod_daily.bmus))

    fig = plt.figure(figsize=[17,12])

    gs3=gridspec.GridSpec(7,7, hspace=0.01, wspace=0.01)
    for b in range(num_clusters):
        s=np.where(sp_mod_daily.bmus==b)[0]
        prob=len(s)/len(sp_mod_daily.bmus)
        ax2=fig.add_subplot(gs3[b],projection='polar')
        z1=4*np.sqrt(np.mean(sp_mod_daily.h.isel(time=s).values,axis=0)) #z1 plots points
        z=np.sum(sp_mod_daily.is_h.isel(time=s).values,axis=0) / len(s) # z is the background
        [p2,p_z1]=plot_spectrum_hs(ax2,np.deg2rad(sp_mod_daily.dir_bins.values), sp_mod_daily.t_bins.values,z, cmap='Greys', vmin=0, vmax=1,alpha_bk=0.4,
                                z1=z1, vmin_z1=0, vmax_z1=vmax_z1, cmap_z1='CMRmap_r', size_point=25,
                                remove_axis=1, prob=prob, prob_max=prob_max, ylim=np.nanmax(sp_mod_daily.t_bins))  
    gs3.tight_layout(fig, rect=[[], [], 0.78, []])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p2,cax=ax0, extend='both')
    ax0.set_ylabel('Cell Probability', fontsize=14)
    gs4.tight_layout(fig, rect=[0.77, 0.1, 0.85, 0.9])

    gs5=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs5[0])
    plt.colorbar(p_z1,cax=ax0, extend='both')
    ax0.set_ylabel('Hs (m)', fontsize=14)
    gs5.tight_layout(fig, rect=[0.85, 0.1, 0.93, 0.9])

    gs6=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs6[0])
    norm = Normalize(vmin=0, vmax=prob_max)        
    cmap = cm.get_cmap('Blues')
    cb1 = mpl.colorbar.ColorbarBase(ax0, cmap=cmap, norm=norm, orientation='vertical', extend='both')
    cb1.set_label('Cluster Probability', fontsize=14)
    gs6.tight_layout(fig, rect=[0.93, 0.1, 1.02, 0.9])
    
def plot_specs_inside_cluster(sp_mod_daily, cluster, gr=6, vmax=0.7):
    
    fig = plt.figure(figsize=[14,12])

    s=np.where(sp_mod_daily.bmus==cluster)[0]

    gs3=gridspec.GridSpec(gr,gr,hspace=0.01, wspace=0.01)
    hs_m=[]

    if len(s)<gr**2:
        lim=len(s)
    else:
        lim=gr**2

    for b in range(np.int(lim)):

        ax2=fig.add_subplot(gs3[b],projection='polar')
        z=sp_mod_daily.h.isel(time=s[b]).values
        z=4*np.sqrt(z)
        hs_m.append(np.sqrt(np.sum(z**2)))
        p2=plot_spectrum_hs(ax2,np.deg2rad(sp_mod_daily.dir_bins.values), sp_mod_daily.t_bins.values,z,vmax=vmax,ylim=np.nanmax(sp_mod_daily.t_bins), remove_axis=1, cmap='magma_r') 

        gs3.tight_layout(fig, rect=[[], [], 0.9, []])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p2,cax=ax0)
    ax0.set_ylabel('Hs (m)', fontsize=14)
    gs4.tight_layout(fig, rect=[0.9, 0.1, 0.99, 0.9])
    
    
def plot_spectrum_hs(ax,x,y,z,z1=[], vmin=0, vmax=0.6,  vmin_z1=0, vmax_z1=0.3, ylim=0.49, size_point=5,point_edge_color=None, alpha_bk=1, cmap='inferno', cmap_z1= 'inferno',remove_axis=0, prob=None, prob_max=0.06, lw=8):
    
    xx,yy=np.meshgrid(x,y)
    
    x=np.append(x,x[0]); 
    y=np.append(0,y); 
    #if vmax<np.nanmax(z):
    #    vmax=1.1*np.nanmax(z)
        
    if cmap=='RdBu_r':
        norm = mcolors.TwoSlopeNorm(0,vmin, vmax)
        p1=ax.pcolormesh(x,y,z, vmin=vmin,vmax=vmax, cmap=plt.cm.RdBu_r, norm=norm, alpha=alpha_bk, linewidths=0.0000001)
    else:
        p1=ax.pcolormesh(x,y,z, vmin=vmin, vmax=vmax,shading='flat', cmap=cmap, alpha=alpha_bk, linewidths=0.00000001)   
    
    if len(z1):
        dx=(x[1]-x[0])/2
        dy=(y[1]-y[0])/2
        p_z1=ax.scatter(xx-dx,yy-dy,size_point, z1, edgecolors=point_edge_color,vmin=vmin_z1,vmax=vmax_z1, cmap=cmap_z1)   
        
    ax.set_theta_zero_location('N', offset=0)
    ax.set_theta_direction(-1)
    ax.set_ylim(0,ylim)
    
    if prob:
        norm = Normalize(vmin=0, vmax=prob_max)        
        cmap = cm.get_cmap('Blues')
        ax.spines['polar'].set_color(cmap(norm(prob)))
        ax.spines['polar'].set_linewidth(lw)
        
    if remove_axis:
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        ax.tick_params(axis='y', colors='plum',labelsize=14,grid_linestyle=':',grid_alpha=0.75,grid_color='plum')
        ax.tick_params(axis='x', colors='purple',labelsize=14,pad=5,grid_linestyle=':',grid_alpha=0.75,grid_color='plum')
    ax.grid(color='plum', linestyle='--', linewidth=0.7,alpha=0.2)
    
    if len(z1):
        return p1, p_z1
    else:
        return p1
    
def plot_forecast_day(slp_d, extent = (50, 290, -67, 65), prob_min=0, prob_max=0.15, s_min=970, s_max=1030):
    
    fig = plt.figure(figsize=[27,7])

    gs3=gridspec.GridSpec(1,1, hspace=0.01, wspace=0.01)

    ax=fig.add_subplot(gs3[0],projection = ccrs.PlateCarree(central_longitude=180))
    ax.set_extent(extent,crs = ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, edgecolor='black',color='grey',zorder=5)
    ax.gridlines()

    vmin, vmax =s_min, s_max
    norm = mcolors.TwoSlopeNorm(1014, vmin=vmin, vmax=vmax)
    z=slp_d.slp_comp.values/100
    z[mask]=np.nan
    p1=ax.pcolormesh(slp_d.longitude.values, slp_d.latitude.values, z,
                      transform=ccrs.PlateCarree(),zorder=2, cmap=plt.cm.RdBu_r, vmin=vmin, vmax=vmax, norm=norm)

    gs3.tight_layout(fig, rect=[[], [], [0.93], []])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p1,cax=ax0, extend='both')
    ax0.set_ylabel('SLP (mbar)', fontsize=14)
    gs4.tight_layout(fig, rect=[0.74, 0.3, 0.79, 0.96])


    gs5=gridspec.GridSpec(1,1, hspace=0.01, wspace=0.01)
    ax2=fig.add_subplot(gs5[0])
    p1=ax2.pcolor(np.flipud((slp_d.prob_spec)), vmin=prob_min, vmax=prob_max, cmap='rainbow')

    gs5.tight_layout(fig, rect=[0.8, 0.25, 1, 1])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p1,cax=ax0, extend='both')
    ax0.set_ylabel('SLP (mbar) - ANOMALIES', fontsize=14)
    gs4.tight_layout(fig, rect=[1, 0.3, 1.05, 0.96])
    
    
def plot_forecast_day_spec(slp_d,spec_kma,mask=[], extent = (50, 290, -67, 65), extent_zoom=[], prob_min=0, prob_max=0.15, s_min=970, s_max=1030, lw=4, pt=[], title=[]):
    
    fig = plt.figure(figsize=[29,10])

    gs3=gridspec.GridSpec(1,1, hspace=0.01, wspace=0.01)

    ax=fig.add_subplot(gs3[0],projection = ccrs.PlateCarree(central_longitude=180))
    ax.set_extent(extent,crs = ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, edgecolor='black',color='grey',zorder=5)
    ax.gridlines()

    vmin, vmax =s_min, s_max
    norm = mcolors.TwoSlopeNorm(1014, vmin=vmin, vmax=vmax)
    z=slp_d.slp_comp.values/100
    if mask:
        z[mask]=np.nan
    p1=ax.pcolormesh(slp_d.longitude.values, slp_d.latitude.values, z,
                      transform=ccrs.PlateCarree(),zorder=2, cmap=plt.cm.RdBu_r, vmin=vmin, vmax=vmax, norm=norm)
    if title:
        ax.text(extent[0]+5, extent[3]-10, title, fontsize=15, fontweight='bold', transform=ccrs.PlateCarree(), zorder=10, color='white')
    
    if pt:
        ax.plot(pt[0],pt[1],'s', markersize=10, color='darkmagenta',transform=ccrs.PlateCarree(), zorder=10)

    gs3.tight_layout(fig, rect=[0, 0.35, 0.6, 1])

    if extent_zoom:
        
        gsz=gridspec.GridSpec(1,1, hspace=0.01, wspace=0.01)

        ax=fig.add_subplot(gsz[0],projection = ccrs.PlateCarree(central_longitude=180))
        ax.set_extent(extent_zoom,crs = ccrs.PlateCarree())
        ax.add_feature(cartopy.feature.LAND, edgecolor='black',color='grey',zorder=5)
        ax.gridlines()

        vmin, vmax =s_min, s_max
        norm = mcolors.TwoSlopeNorm(1014, vmin=vmin, vmax=vmax)
        z=slp_d.slp_comp.values/100
        if mask:
            z[mask]=np.nan
        p1=ax.pcolormesh(slp_d.longitude.values, slp_d.latitude.values, z,
                          transform=ccrs.PlateCarree(),zorder=2, cmap=plt.cm.RdBu_r, vmin=vmin, vmax=vmax, norm=norm)
        if pt:
            ax.plot(pt[0],pt[1],'s', markersize=10, color='darkmagenta',transform=ccrs.PlateCarree(), zorder=10)
        
        plt.colorbar(p1, orientation='horizontal', fraction=0.1, shrink=0.5, pad=0.1).set_label('SLP (mbar)')
        
        gsz.tight_layout(fig, rect=[0, -0.3, 0.6, 0.55])

    num_clusters=len(spec_kma.clusters)

    gs3=gridspec.GridSpec(np.int(np.sqrt(num_clusters)),np.int(np.sqrt(num_clusters)),hspace=0.01, wspace=0.01)
    for b in range(num_clusters):
        s=np.where(spec_kma.bmus==b)[0]
        #prob=len(s)/len(spec_kma.bmus)
        prob=slp_d.prob_spec.values.reshape(-1)[b]
        if prob>0.01:
            ax2=fig.add_subplot(gs3[b],projection='polar')
            z=spec_kma.Mean_h[:,:,b]
            p2=plot_spectrum_hs(ax2,np.deg2rad(spec_kma.dir_bins.values), spec_kma.t_bins.values, z, vmin=0, vmax=1.2, ylim=np.nanmax(spec_kma.t_bins), remove_axis=1,
                                cmap='CMRmap_r', prob=prob, prob_max=prob_max, lw=lw) 
    gs3.tight_layout(fig, rect=[0.55, [], 0.92, []])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p2,cax=ax0)
    ax0.set_ylabel('Hs (m)', fontsize=14)
    gs4.tight_layout(fig, rect=[0.92, 0.1, 0.96, 0.9])

    gs5=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs5[0])
    norm = Normalize(vmin=0, vmax=prob_max)        
    cmap = cm.get_cmap('Blues')
    cb1 = mpl.colorbar.ColorbarBase(ax0, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_label('Probability', fontsize=14)
    gs5.tight_layout(fig, rect=[0.955, 0.1, 1, 0.9])
    
    
def plot_validation_spec(slp_d,spec_kma,mask=[], extent = (50, 290, -67, 65), real_spec=[], prob_min=0, prob_max=0.15, s_min=970, s_max=1030, lw=4, pt=[], title=[], min_prob_clust=0.01):
    
    fig = plt.figure(figsize=[29,10])

    gs3=gridspec.GridSpec(1,1, hspace=0.01, wspace=0.01)

    ax=fig.add_subplot(gs3[0],projection = ccrs.PlateCarree(central_longitude=180))
    ax.set_extent(extent,crs = ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, edgecolor='black',color='grey',zorder=5)
    ax.gridlines()

    vmin, vmax =s_min, s_max
    norm = mcolors.TwoSlopeNorm(1014, vmin=vmin, vmax=vmax)
    z=slp_d.slp_comp.values/100
    if mask:
        z[mask]=np.nan
    p1=ax.pcolormesh(slp_d.longitude.values, slp_d.latitude.values, z,
                      transform=ccrs.PlateCarree(),zorder=2, cmap=plt.cm.RdBu_r, vmin=vmin, vmax=vmax, norm=norm)
    plt.colorbar(p1, orientation='horizontal', fraction=0.1, shrink=0.5, pad=0.1).set_label('SLP (mbar)')
    
    ax.text(extent[0]+5, extent[3]-10, title, fontsize=15, fontweight='bold', transform=ccrs.PlateCarree(), zorder=10, color='white')
    
    if pt:
        ax.plot(pt[0],pt[1],'s', markersize=10, color='darkmagenta',transform=ccrs.PlateCarree(), zorder=10)

    gs3.tight_layout(fig, rect=[0, 0.25, 0.6, 1])

    if len(real_spec):
        
        gsz=gridspec.GridSpec(1,1, hspace=0.01, wspace=0.01)
        
        prob=1
        ax2=fig.add_subplot(gsz[0],projection='polar')
        z=real_spec

        p2=plot_spectrum_hs(ax2,np.deg2rad(spec_kma.dir_bins.values), spec_kma.t_bins.values,z, cmap='CMRmap_r', vmin=0, vmax=1.5, alpha_bk=1,
                                   remove_axis=1, prob=prob, prob_max=prob_max, ylim=np.nanmax(spec_kma.t_bins))  
        
        gsz.tight_layout(fig, rect=[0.3, 0, 0.5, 0.45])

    num_clusters=len(spec_kma.clusters)

    gs3=gridspec.GridSpec(np.int(np.sqrt(num_clusters)),np.int(np.sqrt(num_clusters)),hspace=0.01, wspace=0.01)
    for b in range(num_clusters):
        s=np.where(spec_kma.bmus==b)[0]
        #prob=len(s)/len(spec_kma.bmus)
        prob=slp_d.prob_spec.values.reshape(-1)[b]
        if prob>min_prob_clust:
            ax2=fig.add_subplot(gs3[b],projection='polar')
            z=spec_kma.Mean_h[:,:,b]
            p2=plot_spectrum_hs(ax2,np.deg2rad(spec_kma.dir_bins.values), spec_kma.t_bins.values, z, vmin=0, vmax=1.2, ylim=np.nanmax(spec_kma.t_bins), remove_axis=1,
                                cmap='CMRmap_r', prob=prob, prob_max=prob_max, lw=4) 
    gs3.tight_layout(fig, rect=[0.55, [], 0.92, []])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p2,cax=ax0)
    ax0.set_ylabel('Hs (m)', fontsize=14)
    gs4.tight_layout(fig, rect=[0.92, 0.1, 0.96, 0.9])

    gs5=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs5[0])
    norm = Normalize(vmin=0, vmax=prob_max)        
    cmap = cm.get_cmap('Blues')
    cb1 = mpl.colorbar.ColorbarBase(ax0, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_label('Probability', fontsize=14)
    gs5.tight_layout(fig, rect=[0.955, 0.1, 1, 0.9])


def plot_monthly_mean(slp_dy_m,spec_m,list_months, spec_kma, mean_spec=[], num_clus_pick=49, perct_pick=[], cmap='CMRmap_r', vmin=0, vmax=1.2, alpha_bk=1):

    for month in list_months:

        slp_m=slp_dy_m.sel(time=month).isel(time=0)
        
        if len(mean_spec):
            z=4*np.sqrt(spec_m.sel(time=month).isel(time=0).h.values)- mean_spec
        else:
            z=4*np.sqrt(spec_m.sel(time=month).isel(time=0).h.values)

        fig = plt.figure(figsize=[18,8])

        gs3=gridspec.GridSpec(1,2,hspace=0.1, wspace=0.1)
        ax2=fig.add_subplot(gs3[0],projection='polar')
        ax2.set_ylabel(month, labelpad=40, fontsize=20)
        ax2.set_title('Hindcast', fontsize=20, pad=40)
        p2=plot_spectrum_hs(ax2,np.deg2rad(spec_m.dir_bins.values), spec_m.t_bins.values,z, cmap=cmap, vmin=vmin, vmax=vmax,alpha_bk=alpha_bk,
                            remove_axis=1, prob=1, prob_max=0.2, ylim=np.nanmax(spec_m.t_bins))

        prob=slp_m.prob_spec.values.reshape(-1)
        if num_clus_pick:
            index=np.argsort(-prob)[:num_clus_pick]
            prob=prob[index]
            prob=prob/np.sum(prob)
        elif perct_pick:
            index=np.argsort(-prob)
            prob_s_s=np.cumsum(prob[index])
            index=index[:np.where(prob_s_s>perct_pick)[0][0]]
            prob=prob[index]
            prob=prob/np.sum(prob)

        for b in range(len(index)):

            s=np.where(spec_kma.bmus.values==index[b])[0]
            pr=prob[b]
            h=spec_kma.Mean_h.values[:,:,index[b]]*pr
            if b==0:
                H_c=h
            else:
                H_c=H_c+h

        ax2=fig.add_subplot(gs3[1],projection='polar')
        ax2.set_title('Model', fontsize=20, pad=40)
        if len(mean_spec):
            z1=H_c- mean_spec
        else:
            z1=H_c
            
        p3=plot_spectrum_hs(ax2,np.deg2rad(spec_m.dir_bins.values), spec_m.t_bins.values,z1, cmap=cmap, vmin=vmin, vmax=vmax,alpha_bk=alpha_bk,
                            remove_axis=1, prob=1, prob_max=0.2, ylim=np.nanmax(spec_m.t_bins))  

        gs3.tight_layout(fig, rect=[0.0, 0.1, 0.85, 0.9])

        gs5=gridspec.GridSpec(1,1)
        ax0=fig.add_subplot(gs5[0])
        plt.colorbar(p3,cax=ax0, extend='both')
        ax0.set_ylabel('Hs (m)', fontsize=14)
        gs5.tight_layout(fig, rect=[0.85, 0.2, 0.92, 0.8])

def plot_daily_mean(slp_dy_m,spec_m,list_months, spec_kma, mean_spec=[], num_clus_pick=49, perct_pick=[], cmap='CMRmap_r', vmin=0, vmax=1.2, alpha_bk=1):

    for month in list_months:

        slp_m=slp_dy_m.sel(time=month)
        
        if len(mean_spec):
            z=4*np.sqrt(spec_m.sel(time=month).h.values)- mean_spec
        else:
            z=4*np.sqrt(spec_m.sel(time=month).h.values)

        fig = plt.figure(figsize=[16,6])

        gs3=gridspec.GridSpec(1,2,hspace=0.1, wspace=0.1)
        ax2=fig.add_subplot(gs3[0],projection='polar')
        ax2.set_ylabel(month, labelpad=40, fontsize=20)
        ax2.set_title('Hindcast', fontsize=20, pad=40)
        p2=plot_spectrum_hs(ax2,np.deg2rad(spec_m.dir_bins.values), spec_m.t_bins.values,z, cmap=cmap, vmin=vmin, vmax=vmax,alpha_bk=alpha_bk,
                            remove_axis=1, prob=1, prob_max=0.2, ylim=np.nanmax(spec_m.t_bins))

        prob=slp_m.prob_spec.values.reshape(-1)
        if num_clus_pick:
            index=np.argsort(-prob)[:num_clus_pick]
            prob=prob[index]
            prob=prob/np.sum(prob)
        elif perct_pick:
            index=np.argsort(-prob)
            prob_s_s=np.cumsum(prob[index])
            index=index[:np.where(prob_s_s>perct_pick)[0][0]]
            prob=prob[index]
            prob=prob/np.sum(prob)

        for b in range(len(index)):

            s=np.where(spec_kma.bmus.values==index[b])[0]
            pr=prob[b]
            h=spec_kma.Mean_h.values[:,:,index[b]]*pr
            if b==0:
                H_c=h
            else:
                H_c=H_c+h

        ax2=fig.add_subplot(gs3[1],projection='polar')
        ax2.set_title('Model', fontsize=20, pad=40)
        if len(mean_spec):
            z1=H_c- mean_spec
        else:
            z1=H_c
            
        p3=plot_spectrum_hs(ax2,np.deg2rad(spec_m.dir_bins.values), spec_m.t_bins.values,z1, cmap=cmap, vmin=vmin, vmax=vmax,alpha_bk=alpha_bk,
                            remove_axis=1, prob=1, prob_max=0.2, ylim=np.nanmax(spec_m.t_bins))  

        gs3.tight_layout(fig, rect=[0.0, 0.1, 0.85, 0.9])

        gs5=gridspec.GridSpec(1,1)
        ax0=fig.add_subplot(gs5[0])
        plt.colorbar(p3,cax=ax0, extend='both')
        ax0.set_ylabel('Hs (m)', fontsize=14)
        gs5.tight_layout(fig, rect=[0.85, 0.2, 0.92, 0.8])

def plot_spec_n34_mjo_month(sp_mod_daily, n34_sel=[], mjo_sel=[], month=[], an_range=0.1, mean_month=[], figsize=[22,7], vmax_z1=0.8, vmax_an=0.07):

    fig = plt.figure(figsize=figsize)
    
    gr=6
    gs3=gridspec.GridSpec(1,2,hspace=0.1, wspace=0.1)
    
    if n34_sel and mjo_sel and month:
        label= 'n34= ' + str(n34_sel)+ ' |  MJO= ' + str(mjo_sel) + ' |  Month= ' + str(month)
        s=np.where((sp_mod_daily.n34==n34_sel-1) & (sp_mod_daily.mjo==mjo_sel) & (sp_mod_daily.time.dt.month.values==month))[0]
    elif n34_sel and mjo_sel:
        label= 'n34= ' + str(n34_sel)+ ' |  MJO= ' + str(mjo_sel) 
        s=np.where((sp_mod_daily.n34==n34_sel-1) & (sp_mod_daily.mjo==mjo_sel))[0]
    elif n34_sel and month:
        label= 'n34= ' + str(n34_sel)+ ' |  Month= ' + str(month)
        s=np.where((sp_mod_daily.n34==n34_sel-1) & (sp_mod_daily.time.dt.month.values==month))[0]
    elif mjo_sel and month:
        label= ' |  MJO= ' + str(mjo_sel) + ' |  Month= ' + str(month)
        s=np.where((sp_mod_daily.mjo==mjo_sel) & (sp_mod_daily.time.dt.month.values==month))[0]
    elif mjo_sel:
        label=  ' MJO= ' + str(mjo_sel) 
        s=np.where(sp_mod_daily.mjo==mjo_sel)[0]
    elif n34_sel:
        label= 'n34= ' + str(n34_sel)
        s=np.where((sp_mod_daily.n34==n34_sel-1))[0]
    elif month:
        label=' Month= ' + str(month)
        s=np.where((sp_mod_daily.time.dt.month.values==month))[0]
    else:
        label='ALL DATA'
        s=np.arange(len(sp_mod_daily.time))
    
    print('Number of data: ' + str(len(s)))
    prob=len(s)/len(sp_mod_daily.time)
    ax2=fig.add_subplot(gs3[0],projection='polar')
    ax2.set_ylabel(label, labelpad=40, fontsize=20)

    z1=np.mean(sp_mod_daily.h.isel(time=s).values,axis=0) #z1 plots points
    z1=4*np.sqrt(z1)
    z=np.sum(sp_mod_daily.is_h.isel(time=s).values,axis=0) / len(s) # z is the background

    [p2,p_z1]=plot_spectrum_hs(ax2,np.deg2rad(sp_mod_daily.dir_bins.values), sp_mod_daily.t_bins.values,z, cmap='Greys', vmin=0, vmax=1, alpha_bk=0.4,
                            z1=z1, vmin_z1=0, vmax_z1=vmax_z1, cmap_z1='CMRmap_r', size_point=100,
                            remove_axis=1, prob=prob, prob_max=0.2, ylim=np.nanmax(sp_mod_daily.t_bins))  

    ax2=fig.add_subplot(gs3[1],projection='polar')   
    
    if mean_month:
        s1=np.where((sp_mod_daily.time.dt.month.values==month))[0]
        mean_spec=4*np.sqrt(np.mean(sp_mod_daily.isel(time=s1).h.values,axis=0))
        mean_prob=np.sum(sp_mod_daily.isel(time=s1).is_h.values,axis=0) / len(sp_mod_daily.isel(time=s1).time)
    else:
        mean_spec=4*np.sqrt(np.mean(sp_mod_daily.h.values,axis=0))
        mean_prob=np.sum(sp_mod_daily.is_h.values,axis=0) / len(sp_mod_daily.time)
    z1 = z1 - mean_spec #z1 plots points 
    z = z - mean_prob
    
    [p3,p_z1_an]=plot_spectrum_hs(ax2,np.deg2rad(sp_mod_daily.dir_bins.values), sp_mod_daily.t_bins.values,z,point_edge_color='Grey', cmap='RdBu_r', vmin=-vmax_an, vmax=vmax_an ,alpha_bk=0.3,
                            z1=z1, vmin_z1=-an_range, vmax_z1=an_range, cmap_z1='RdBu_r', size_point=90,
                            remove_axis=1, prob=prob, prob_max=0.2, ylim=np.nanmax(sp_mod_daily.t_bins))  

    gs3.tight_layout(fig, rect=[0, [], 0.8, []])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p2,cax=ax0, extend='both')
    ax0.set_ylabel('Cell Probability', fontsize=14)
    gs4.tight_layout(fig, rect=[0.755, 0.1, 0.82, 0.9])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p3,cax=ax0, extend='both')
    ax0.set_ylabel('Cell Probability - Anomaly', fontsize=14)
    gs4.tight_layout(fig, rect=[0.82, 0.1, 0.89, 0.9])

    gs5=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs5[0])
    plt.colorbar(p_z1,cax=ax0, extend='both')
    ax0.set_ylabel('Hs (m)', fontsize=14)
    gs5.tight_layout(fig, rect=[0.89, 0.1, 0.95, 0.9])

    gs6=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs6[0])
    plt.colorbar(p_z1_an,cax=ax0, extend='both')
    ax0.set_ylabel('Hs (m) - Amomaly', fontsize=14)
    gs6.tight_layout(fig, rect=[0.95, 0.1, 1.02, 0.9])

def spec_n34_mjo_month_ax(sp_mod_daily, n34_sel=[], mjo_sel=[], month=[], an_range=0.1, ax2=[], anomaly=[], mean_month=[], size_point=18, lw=3, vmax_z1=0.7):

   
    if n34_sel and mjo_sel and month:
        label= 'n34= ' + str(n34_sel)+ ' |  MJO= ' + str(mjo_sel) + ' |  Month= ' + str(month)
        s=np.where((sp_mod_daily.n34==n34_sel-1) & (sp_mod_daily.mjo==mjo_sel) & (sp_mod_daily.time.dt.month.values==month))[0]
    elif n34_sel and mjo_sel:
        label= 'n34= ' + str(n34_sel)+ ' |  MJO= ' + str(mjo_sel) 
        s=np.where((sp_mod_daily.n34==n34_sel-1) & (sp_mod_daily.mjo==mjo_sel))[0]
    elif n34_sel and month:
        label= 'n34= ' + str(n34_sel)+ ' |  Month= ' + str(month)
        s=np.where((sp_mod_daily.n34==n34_sel-1) & (sp_mod_daily.time.dt.month.values==month))[0]
    elif mjo_sel and month:
        label= ' |  MJO= ' + str(mjo_sel) + ' |  Month= ' + str(month)
        s=np.where((sp_mod_daily.mjo==mjo_sel) & (sp_mod_daily.time.dt.month.values==month))[0]
    elif mjo_sel:
        label=  ' MJO= ' + str(mjo_sel) 
        s=np.where(sp_mod_daily.mjo==mjo_sel)[0]
    elif n34_sel:
        label= 'n34= ' + str(n34_sel)
        s=np.where((sp_mod_daily.n34==n34_sel-1))[0]
    elif month:
        label=' Month= ' + str(month)
        s=np.where((sp_mod_daily.time.dt.month.values==month))[0]
    else:
        label='ALL DATA'
        s=np.arange(len(sp_mod_daily.time))
    
    prob=len(s)/len(sp_mod_daily.time)
    
    z1=4*np.sqrt(np.mean(sp_mod_daily.h.isel(time=s).values,axis=0)) #z1 plots points
    z=np.sum(sp_mod_daily.is_h.isel(time=s).values,axis=0) / len(s) # z is the background

    
    if anomaly:
        if mean_month:
            s1=np.where((sp_mod_daily.time.dt.month.values==month))[0]
            mean_spec=4*np.sqrt(np.mean(sp_mod_daily.isel(time=s1).h.values,axis=0))
            mean_prob=np.sum(sp_mod_daily.isel(time=s1).is_h.values,axis=0) / len(sp_mod_daily.isel(time=s1).time)
        else:
            mean_spec=4*np.sqrt(np.mean(sp_mod_daily.h.values,axis=0))
            mean_prob=np.sum(sp_mod_daily.is_h.values,axis=0) / len(sp_mod_daily.time)
        z1 = z1 - mean_spec #z1 plots points
        z = z - mean_prob
        [p2,p_z1]=plot_spectrum_hs(ax2,np.deg2rad(sp_mod_daily.dir_bins.values), sp_mod_daily.t_bins.values,z,point_edge_color='Grey', cmap='RdBu_r', 
                            vmin=-0.04, vmax=0.04 ,alpha_bk=0.3, lw=lw,
                            z1=z1, vmin_z1=-an_range, vmax_z1=an_range, cmap_z1='RdBu_r', size_point=size_point,
                            remove_axis=1, prob=prob, prob_max=0.02, ylim=np.nanmax(sp_mod_daily.t_bins))
    else:
        [p2,p_z1]=plot_spectrum_hs(ax2,np.deg2rad(sp_mod_daily.dir_bins.values), sp_mod_daily.t_bins.values,z, cmap='Greys', vmin=0, vmax=1, alpha_bk=0.4,
                            z1=z1, vmin_z1=0, vmax_z1=vmax_z1, cmap_z1='CMRmap_r', size_point=size_point, lw=lw,
                            remove_axis=1, prob=prob, prob_max=0.02, ylim=np.nanmax(sp_mod_daily.t_bins))
        
    
    
    return p2, p_z1


def plot_spec_n34(sp_mod_daily, annomaly=[], vmax_z1=0.7):
    
    fig = plt.figure(figsize=[18,4])

    gr=3
    name=['El Niño', 'Neutral', 'La Niña']
    gs3=gridspec.GridSpec(1,3,hspace=0.01, wspace=0.01)
    for b in range(gr):

        s=np.where(sp_mod_daily.n34==b)[0]
        prob=len(s)/len(sp_mod_daily.time)
        ax2=fig.add_subplot(gs3[b],projection='polar')
        ax2.set_title(name[b], fontsize=16)
        if annomaly:
            mean_spec=4*np.sqrt(np.mean(sp_mod_daily.h.values,axis=0)) #z1 plots points
            mean_prob=np.sum(sp_mod_daily.is_h.values,axis=0) / len(sp_mod_daily.time) # z is the background
            z1=4*np.sqrt(np.mean(sp_mod_daily.h.isel(time=s).values,axis=0))
            z1=z1 - mean_spec #z1 plots points
            z=np.sum(sp_mod_daily.is_h.isel(time=s).values,axis=0) / len(s) - mean_prob # z is the background
            prob_max=0.3
            [p2,p_z1]=plot_spectrum_hs(ax2,np.deg2rad(sp_mod_daily.dir_bins.values), sp_mod_daily.t_bins.values,z,point_edge_color='Grey', cmap='RdBu_r', vmin=-0.04, vmax=0.04,alpha_bk=0.4,
                                    z1=z1, vmin_z1=-vmax_z1, vmax_z1=vmax_z1, cmap_z1='RdBu_r', size_point=50, lw=5,
                                    remove_axis=1, prob=prob, prob_max=prob_max, ylim=np.nanmax(sp_mod_daily.t_bins))  
        else:
            z1=4*np.sqrt(np.mean(sp_mod_daily.h.isel(time=s).values,axis=0)) #z1 plots points
            z=np.sum(sp_mod_daily.is_h.isel(time=s).values,axis=0) / len(s) # z is the background
            prob_max=0.2
            [p2,p_z1]=plot_spectrum_hs(ax2,np.deg2rad(sp_mod_daily.dir_bins.values), sp_mod_daily.t_bins.values,z, cmap='Greys', vmin=0, vmax=1,alpha_bk=0.4,
                                    z1=z1, vmin_z1=0, vmax_z1=vmax_z1, cmap_z1='CMRmap_r', size_point=40,lw=5,
                                    remove_axis=1, prob=prob, prob_max=prob_max, ylim=np.nanmax(sp_mod_daily.t_bins))  
    gs3.tight_layout(fig, rect=[[], [], 0.78, []])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p2,cax=ax0, extend='both')
    ax0.set_ylabel('Cell Probability', fontsize=14)
    gs4.tight_layout(fig, rect=[0.77, 0.1, 0.85, 0.9])

    gs5=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs5[0])
    plt.colorbar(p_z1,cax=ax0, extend='both')
    ax0.set_ylabel('Hs (m)', fontsize=14)
    gs5.tight_layout(fig, rect=[0.85, 0.1, 0.93, 0.9])

    gs6=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs6[0])
    norm = Normalize(vmin=0, vmax=prob_max)        
    cmap = cm.get_cmap('Blues')
    cb1 = mpl.colorbar.ColorbarBase(ax0, cmap=cmap, norm=norm, orientation='vertical', extend='both')
    cb1.set_label('Cluster Probability', fontsize=14)
    gs6.tight_layout(fig, rect=[0.93, 0.1, 1.02, 0.9])
    
    
def plot_spec_mjo(sp_mod_daily, annomaly=[], vmax_z1=0.7):
    
    fig = plt.figure(figsize=[19,7])

    gr=8
    gs3=gridspec.GridSpec(2,4,hspace=0.01, wspace=0.01)
    hs_m=[]

    for b in range(gr):

        s=np.where(sp_mod_daily.mjo==b+1)[0]
        prob=len(s)/len(sp_mod_daily.bmus)
        ax2=fig.add_subplot(gs3[b],projection='polar')
        ax2.set_title('MJO' + str(b+1), fontsize=16)
        if annomaly:
            
            mean_spec=4*np.sqrt(np.mean(sp_mod_daily.h.values,axis=0)) #z1 plots points
            mean_prob=np.sum(sp_mod_daily.is_h.values,axis=0) / len(sp_mod_daily.time) # z is the background
            
            z1=4*np.sqrt(np.mean(sp_mod_daily.h.isel(time=s).values,axis=0)) - mean_spec #z1 plots points
            z=np.sum(sp_mod_daily.is_h.isel(time=s).values,axis=0) / len(s)  - mean_prob# z is the background
            
            prob_max=0.15
    
            [p2,p_z1]=plot_spectrum_hs(ax2,np.deg2rad(sp_mod_daily.dir_bins.values), sp_mod_daily.t_bins.values,z, point_edge_color='Grey', cmap='RdBu_r', vmin=-0.04, vmax=0.04,
                            alpha_bk=0.4, z1=z1, vmin_z1=-vmax_z1, vmax_z1=vmax_z1, cmap_z1='RdBu_r', size_point=45, lw=5,
                            remove_axis=1, prob=prob, prob_max=prob_max, ylim=np.nanmax(sp_mod_daily.t_bins))
        else:
            z1=4*np.sqrt(np.mean(sp_mod_daily.h.isel(time=s).values,axis=0)) #z1 plots points
            z=np.sum(sp_mod_daily.is_h.isel(time=s).values,axis=0) / len(s) # z is the background
            prob_max=0.2
            [p2,p_z1]=plot_spectrum_hs(ax2,np.deg2rad(sp_mod_daily.dir_bins.values), sp_mod_daily.t_bins.values,z, cmap='Greys', vmin=0, vmax=1,alpha_bk=0.4,
                                z1=z1, vmin_z1=0, vmax_z1=vmax_z1, cmap_z1='CMRmap_r', size_point=45, lw=5,
                                remove_axis=1, prob=prob, prob_max=prob_max, ylim=np.nanmax(sp_mod_daily.t_bins))  

    gs3.tight_layout(fig, rect=[[], [], 0.78, []])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p2,cax=ax0, extend='both')
    ax0.set_ylabel('Cell Probability', fontsize=14)
    gs4.tight_layout(fig, rect=[0.77, 0.1, 0.85, 0.9])

    gs5=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs5[0])
    plt.colorbar(p_z1,cax=ax0, extend='both')
    ax0.set_ylabel('Hs (m)', fontsize=14)
    gs5.tight_layout(fig, rect=[0.85, 0.1, 0.93, 0.9])

    gs6=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs6[0])
    norm = Normalize(vmin=0, vmax=prob_max)        
    cmap = cm.get_cmap('Blues')
    cb1 = mpl.colorbar.ColorbarBase(ax0, cmap=cmap, norm=norm, orientation='vertical', extend='both')
    cb1.set_label('Cluster Probability', fontsize=14)
    gs6.tight_layout(fig, rect=[0.93, 0.1, 1.02, 0.9])
    
    
def plot_slp_field(slp, extent=[], time=4804):
    
    fig = plt.figure(figsize=(24, 10))

    ax = plt.axes(projection = ccrs.PlateCarree(central_longitude=180))
    if extent:
        ax.set_extent(extent,crs = ccrs.PlateCarree())
    ax.stock_img()

    # cartopy land feature
    land_10m = cartopy.feature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='darkgrey', facecolor='gainsboro',  zorder=5)
    ax.add_feature(land_10m)
    ax.gridlines()

    # scatter data points
    vmin, vmax =96000, 103000

    norm = mcolors.TwoSlopeNorm(101400)
    p1=plt.pcolor(slp.longitude.values, slp.latitude.values, slp.slp.isel(time=time), transform=ccrs.PlateCarree(),zorder=2, cmap=plt.cm.RdBu_r, norm=norm)
    plt.colorbar(p1).set_label('SLP')
    plt.title('SLP ' + str(slp.isel(time=4804).time.values));

def plot_slp_kmeans(slp, extent=(50, 290, -67, 65),figsize=[15.5,14], dynamic=[], annomaly=[]):
    
    fig = plt.figure(figsize=figsize)
    extent = (50, 290, -67, 65)

    mask=np.where(np.isnan(slp.est_mask)==True)
    num_clusters=len(np.unique(slp.bmus))

    gr=np.int(np.sqrt(num_clusters))
    gs3=gridspec.GridSpec(gr,gr, hspace=0.01, wspace=0.01)
    for b in range(num_clusters):

        s=np.where(slp.bmus==b)[0]
        ax2=fig.add_subplot(gs3[b],projection = ccrs.PlateCarree(central_longitude=180))
        ax2.set_extent(extent,crs = ccrs.PlateCarree())

        ax2.add_feature(cartopy.feature.LAND, edgecolor='black',color='grey',zorder=5)
        ax2.gridlines()

        if annomaly:
            vmin, vmax =-20, 10
            norm = mcolors.TwoSlopeNorm(0, vmin=vmin, vmax=vmax)
            mean_slp=np.nanmean(slp.slp.values, axis=0)/100
            z=np.nanmean(slp.slp.isel(time=s).values/100-mean_slp, axis=0)    
        else:
            vmin, vmax =990, 1030
            norm = mcolors.TwoSlopeNorm(1014, vmin=vmin, vmax=vmax)
            z=np.nanmean(slp.slp.isel(time=s).values, axis=0)/100
        z[mask]=np.nan
        p1=ax2.pcolormesh(slp.longitude.values, slp.latitude.values, z, transform=ccrs.PlateCarree(),zorder=2, cmap=plt.cm.RdBu_r, vmin=vmin, vmax=vmax, norm=norm)


    gs3.tight_layout(fig, rect=[[], [], 0.93, []])
    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p1,cax=ax0, extend='both')
    ax0.set_ylabel('SLP (mbar)', fontsize=14)
    gs4.tight_layout(fig, rect=[0.92, 0.15, 0.985, 0.85])
    
def plot_entropy(slp_dy_m, re_all):
    
    plt.figure(figsize=[20,4])
    plt.grid('both', color='plum', linestyle='--')
    plt.plot(slp_dy_m.time.values[:-4], re_all[1:-1], '.-', linewidth=2, markersize=15, color='firebrick')
    plt.xlabel('Time', fontsize=15)
    plt.ylabel('RE', fontsize=15)
    plt.xlim([slp_dy_m.time.values[1],slp_dy_m.time.values[-4]])
    plt.text(slp_dy_m.time.values[-11],1.4,'Mean RE = ' + str(np.round(np.nanmean(re_all[1:]),2)), fontsize=15, color='darkred', fontweight='bold');
    
    
def plot_prob_comparison(sp_mod_daily,slp_dy_m,spec_kma, month):
    
    num_clusters_spec=len(np.unique(sp_mod_daily.bmus))

    sp=sp_mod_daily.sel(time=month[0])
    real_prob=np.full([len(spec_kma.clusters)],0.0)
    for a in range(len(spec_kma.clusters)):
        s=np.where(sp.bmus.values==a)[0]
        real_prob[a]=len(s)/len(sp.time)

    slp_m=slp_dy_m.sel(time=month[0]).isel(time=0)
    model_prob=slp_m.prob_spec.values.reshape(-1)

    real_prob[np.where(real_prob==0)]=0.00000000000000001
    model_prob[np.where(model_prob==0)]=0.00000000000000001
    RE=entropy(real_prob, model_prob)
    print('RE= ' + str(RE))

    prob_real=np.reshape(real_prob,[np.sqrt(num_clusters_spec).astype('int'), np.sqrt(num_clusters_spec).astype('int')])
    prob_model=np.reshape(model_prob,[np.sqrt(num_clusters_spec).astype('int'), np.sqrt(num_clusters_spec).astype('int')])

    fig = plt.figure(figsize=[19,7])
    gs3=gridspec.GridSpec(1,3, hspace=0.05, wspace=0.1)

    ax=fig.add_subplot(gs3[0])
    plt.pcolor(np.flipud(prob_real), cmap='Blues', vmin=0, vmax=0.25)
    ax.set_title('Real', fontsize=15)
    plt.colorbar(orientation='horizontal', pad=0.06, extend='both').set_label('Prob', fontsize=15)

    ax=fig.add_subplot(gs3[1])
    plt.pcolor(np.flipud(prob_model), cmap='Blues', vmin=0, vmax=0.25)
    ax.set_title('Model', fontsize=15)
    plt.colorbar(orientation='horizontal', pad=0.06, extend='both').set_label('Prob', fontsize=15)

    ax=fig.add_subplot(gs3[2])
    plt.pcolor(np.flipud(prob_real - prob_model), cmap='RdBu_r', vmin=-0.25, vmax=0.25)
    ax.set_title('Real - Model', fontsize=15)
    plt.colorbar(orientation='horizontal', pad=0.06, extend='both').set_label('\u0394Prob', fontsize=15)
    
    
def plot_monthly_mean_ax(spec_kma,slp_dy_m, list_months, ax2, mean_spec=[], num_clus_pick=49, perct_pick=[], cmap='CMRmap_r', vmin=0, vmax=1.2, alpha_bk=1):

    for month in list_months:
        
        slp_m=slp_dy_m.sel(time=month).isel(time=0)
        prob=slp_m.prob_spec.values.reshape(-1)
        
        if num_clus_pick:
            index=np.argsort(-prob)[:num_clus_pick]
            prob=prob[index]
            prob=prob/np.sum(prob)
        elif perct_pick:
            index=np.argsort(-prob)
            prob_s_s=np.cumsum(prob[index])
            index=index[:np.where(prob_s_s>perct_pick)[0][0]]
            prob=prob[index]
            prob=prob/np.sum(prob)

        for b in range(len(index)):

            s=np.where(spec_kma.bmus.values==index[b])[0]
            pr=prob[b]
            h=spec_kma.Mean_h.values[:,:,index[b]]*pr
            if b==0:
                H_c=h
            else:
                H_c=H_c+h

        ax2.set_title(month, fontsize=15, pad=10)
        if len(mean_spec):
            z1=H_c- mean_spec
        else:
            z1=H_c
            
        p3=plot_spectrum_hs(ax2,np.deg2rad(spec_kma.dir_bins.values), spec_kma.t_bins.values,z1, cmap=cmap, vmin=vmin, vmax=vmax,alpha_bk=alpha_bk,
                            remove_axis=1, prob=1, prob_max=0.2, ylim=np.nanmax(spec_kma.t_bins), lw=3) 
        
    return p3
       
def Plot_n34_forecast(n34_f, l1, l2):
    plt.figure(figsize=[15,6])
    plt.grid(color='plum')

    for ens in n34_f.ENS:
        plt.plot(n34_f.TIME,n34_f.anom.sel(ENS=ens),'-', linewidth=1.5, markersize=8, alpha=0.5, zorder=2)

    ensemble_f=np.nanmean(n34_f.anom, axis=0)
    plt.plot(n34_f.TIME[1:],ensemble_f[1:],'-', color='navy', linewidth=3, markersize=8, alpha=0.8, zorder=3)

    l1=0.5
    l2=-0.5

    clasif=np.full([len(n34_f.TIME)],1.0)
    s1=np.where(ensemble_f>l1)[0]
    s2=np.where(ensemble_f<l2)[0]
    s3=np.where((ensemble_f>l2) & (ensemble_f<l1))[0]

    clasif[s1]=0
    clasif[s2]=2
    plt.plot(n34_f.TIME[s1],ensemble_f[s1],'.',color='indianred', markersize=25,label='El Niño [1]', zorder=4)
    plt.plot(n34_f.TIME[s3],ensemble_f[s3],'.',color='darkgreen', markersize=25,label='Neutral [2]', zorder=4)
    plt.plot(n34_f.TIME[s2],ensemble_f[s2],'.',color='cornflowerblue', markersize=25,label='La Niña [3]', zorder=4)

    n34_f['ensemble']=('time',ensemble_f) #0:Niño, 1:Neutral, 2:Niña
    n34_f['class_ensemble']=('time',clasif) #0:Niño, 1:Neutral, 2:Niña

    plt.xlabel('Time',fontsize=16)
    plt.ylabel('Niño 3.4',fontsize=16)

    plt.plot([n34_f.TIME.values[0], n34_f.TIME.values[-1]], [l1,l1],'-',color='indianred', linewidth=2, label='El Niño threshold')
    plt.plot([n34_f.TIME.values[0], n34_f.TIME.values[-1]], [0,0],'--',color='grey', linewidth=2)
    plt.plot([n34_f.TIME.values[0], n34_f.TIME.values[-1]], [l2,l2],'-',color='cornflowerblue', linewidth=2, label='La Niña threshold')

    plt.xlim([n34_f.TIME[1], n34_f.TIME[-1]])
    plt.legend(ncol=5, fontsize=13)
    
    return n34_f
