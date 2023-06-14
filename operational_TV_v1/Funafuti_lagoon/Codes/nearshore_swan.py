import warnings
warnings.filterwarnings('ignore')

import sys
import os
import time
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cmocean
from scipy import interpolate
from scipy.stats import gaussian_kde

import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize,ListedColormap, LinearSegmentedColormap
import matplotlib as mpl
from matplotlib import cm
from matplotlib import animation
import pickle as pk

from sklearn.cluster import KMeans
from scipy.stats import circmean

from math import gamma as gm
import wavespectra

import cartopy
import cartopy.crs as ccrs
from cartopy.feature import ShapelyFeature
from cartopy.io.shapereader import Reader


#############################################################
######################   CALCULATIONS  ######################
#############################################################


def interpolate_bathymetry(bati_cut, resolution):
    
    '''
    bati_cut :   original bathymetry
    resolution:  resolution for the interpolation
    '''
    
    xn=np.arange(np.min(bati_cut.x), np.nanmax(bati_cut.x), resolution)
    yn=np.arange(np.min(bati_cut.y), np.nanmax(bati_cut.y), resolution)
    
    f = interpolate.interp2d(bati_cut.x, bati_cut.y, bati_cut.z, kind='linear')
    znew = f(xn,yn)
    bati_interp = xr.Dataset(
        {
            'z': (['y','x'], znew),
        },
        coords={
            'y': yn,
            'x': xn,
        },
    )
    
    return bati_interp

def Reconstruct_forecast(ws, dirs_rec, freqs_rec, gamma_swell=3.3 , gamma_sea=3.3, dspr_sea=30.0, dspr_swell=10.0):
    
    '''
    This function reconstruct the spectral forecast from the different partitions.
    
    ws         :  Forecast partitions
    dirs_rec   :  Spec directions for reconstruction
    freqs_rec  :  Frequencies for reconstruction
    
    '''

    Sum_Snn = np.zeros((len(dirs_rec),len(freqs_rec),len(ws.time)))
    print('Times: ' + str(len(ws.time)))

    for rr in range(len(ws.time)):
        
        #sys.stdout.flush()
        #sys.stdout.write('\rPos: %d\r' %rr)
        #sys.stdout.flush()

        ws1=ws.isel(time=rr)
        hh=np.stack((ws1.hs0.values,ws1.hs1.values, ws1.hs2.values,ws1.hs3.values))
        tt=np.stack((ws1.tp0.values, ws1.tp1.values, ws1.tp2.values,ws1.tp3.values))
        dd=np.stack((ws1.dir0.values,ws1.dir1.values, ws1.dir2.values,ws1.dir3.values))
        gamma=np.array((gamma_sea, gamma_swell, gamma_swell, gamma_swell))
        dspr=np.deg2rad(np.array((dspr_sea,dspr_swell,dspr_swell,dspr_swell)))
        s = (2/(dspr**2))-1

        tt=tt[~np.isnan(hh)]
        dd=dd[~np.isnan(hh)]
        gamma=gamma[~np.isnan(hh)]
        s=s[~np.isnan(hh)]
        hh=hh[~np.isnan(hh)]

        for p in range(len(hh)):                       

            sigma=np.full([len(freqs_rec),1],0.07)
            sigma[np.where(freqs_rec>(1/ws1['tp'+str(p)].values))]=0.09

            Beta = (0.06238/(0.23+0.0336*gamma[p]-0.185*(1.9+gamma[p])**-1))*(1.094-0.01915*np.log(gamma[p]))
            S = Beta * (ws1['hs'+str(p)].values**2) * (ws1['tp'+str(p)].values**-4) * (freqs_rec**-5)*np.exp(-1.25*(ws1['tp'+str(p)].values*freqs_rec)**-4) \
                *gamma[p]**(np.exp((-(ws1['tp'+str(p)].values*freqs_rec-1)**2)/(2*sigma.T**2)))
            D = ((2**(2*s[p]-1))/np.pi)*(gm(s[p]+1)**2/gm(2*s[p]+1))*np.abs(np.cos((np.deg2rad(dirs_rec)-np.deg2rad(ws1['dir'+str(p)].values))/2))**(2*s[p])
            Snn = np.multiply(S.T,D) * (np.pi/180)
            Snn[np.isnan(Snn)] = 0
            Sum_Snn[:,:,rr] = Sum_Snn[:,:,rr] + Snn.T    

    spec_or = xr.Dataset({'efth': (['dir','freq','time'],Sum_Snn)}, coords={'dir': dirs_rec, 'freq': freqs_rec, 'time':ws.time})
    st_or = spec_or.spec.stats(['hs','tp','tm02','dpm','dspr'])
        
    return spec_or, st_or

def process_waves_f(waves):
    
    waves_f = xr.Dataset({'hs0': (['time','lat','lon'],waves.wvhgtsfc.values),
                         'hs1': (['time','lat','lon'],waves.swell_1.values),
                         'hs2': (['time','lat','lon'],waves.swell_2.values),
                         'hs3': (['time','lat','lon'],waves.swell_3.values),
                         'dir0': (['time','lat','lon'],waves.wvdirsfc.values),
                         'dir1': (['time','lat','lon'],waves.swdir_1.values),
                         'dir2': (['time','lat','lon'],waves.swdir_2.values),
                         'dir3': (['time','lat','lon'],waves.swdir_3.values),
                         'tp0': (['time','lat','lon'],1.3*waves.wvpersfc.values),
                         'tp1': (['time','lat','lon'],1.3*waves.swper_1.values),
                         'tp2': (['time','lat','lon'],1.3*waves.swper_2.values),
                         'tp3': (['time','lat','lon'],1.3*waves.swper_3.values),
                         'Wdir': (['time','lat','lon'],waves.wdirsfc.values),
                         'Wspeed': (['time','lat','lon'],waves.windsfc.values),
                          
                   }, 
                  coords={'time': waves.time.values, 
                          'lat': waves.lat.values, 
                          'lon': waves.lon.values})
    return waves_f

def spec_kmeans_clustering_NS(path_main,site, sp_mod_daily, n_pcs, num_clusters, do=[], min_data_frac=15):
    
    if do==1:
        min_data=np.int(len(sp_mod_daily.time)/num_clusters/min_data_frac)
        num, it=0,0

        while np.nanmin(num)<min_data:

            sys.stdout.write('\rIteration: %d' %it)
            sys.stdout.flush()

            kma = KMeans(n_clusters=num_clusters,init='k-means++',verbose=0, n_init=5, tol=0.001) #80
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

        pk.dump(kma, open(os.path.join(path_main, 'Results', site , 'Spec_KMA_' + site + '_NS.pkl'),"wb"))

        print('Number of sims: ' + str(kma.n_iter_))
        _, group_size = np.unique(bmus, return_counts=True)
        print('Minimum number of data: ' + str(np.nanmin(group_size)))
        sp_mod_daily['kma_order']=kma_order
        sp_mod_daily['n_pcs']=n_pcs
        sp_mod_daily['bmus']=('time',bmus)
        sp_mod_daily['centroids']=(('clusters','n_pcs_'),centers)
        sp_mod_daily.to_netcdf(os.path.join(path_main,'Results',site,'Spec_classification_KMA_' + str(num_clusters) + '_NS.nc'))
    else:
        sp_mod_daily=xr.open_dataset(os.path.join(path_main,'Results',site,'Spec_classification_KMA_' + str(num_clusters) + '_NS.nc'))

    return sp_mod_daily

def Calculate_mean_clusters(sp, annomaly=[]):
    
    num_clusters=len(np.unique(sp.bmus))
    mean_z=np.mean(sp.efth.values,axis=0)
    
    h_clust=np.full([len(sp.freq), len(sp.dir), len(sp.kma_order)], np.nan)

    for b in range(num_clusters):
        s=np.where(sp.bmus==b)[0]
        prob=len(s)/len(sp.bmus)
        if annomaly:
            z=np.mean(sp.efth.isel(time=s).values,axis=0) - mean_z
            h_clust[:,:,b]=z
        else:
            z=np.mean(sp.efth.isel(time=s).values,axis=0)
            h_clust[:,:,b]=z
            
    
    return h_clust

def Generate_kp_coefficients(save_kp_path, out_sim, sp_kma):
    
    '''
    This function processes the propagation coefficients for all the grid points within the simulation output (out_sim).
    It takes a long time to run but it only needs to be done once per location
    save_kp_path    :    location to save all the kp coefficients. One file per location
    out_sim         :    results from simulations. Dimensions: cases x lat x lon
    sp_kma          :    frequency and direction coordinates
    
    '''
    
    print('Preprocessing and saving kp coefficients for all the grid points')
    print('Total: ' + str(len(out_sim.lon) * len(out_sim.lat)))
    start=time.time()
    
    cc=0
    ini=time.time()

    for ilon in np.arange(0,len(out_sim.lon),1):
        for ilat in range(len(out_sim.lat)):

            sys.stdout.write('\rCount %d' %cc)
            sys.stdout.flush()
            cc+=1

            kp=np.full([len(out_sim.case), len(sp_kma.freq.values), len(sp_kma.dir.values)],0.0)

            for case in out_sim.case:

                out_sim_pc=out_sim.isel(lon=ilon, lat=ilat, case=case)
                mask=np.where((out_sim_pc.TPsmoo_part.values>0))[0]
                k=out_sim_pc.Hs_part.values[mask]
                f=1/out_sim_pc.TPsmoo_part.values[mask]
                d=out_sim_pc.Dir_part.values[mask]

                for c in range(len(f)):
                    i=np.argmin(np.abs(sp_kma.freq.values-f[c]))
                    j=np.argmin(np.abs(sp_kma.dir.values-d[c]))
                    kp[case, i, j]=k[c] #TODO: Cuidado si caen dos en la misma celda porq solo coge el ultimo

            pk.dump(kp, open(os.path.join(save_kp_path, 'kp_lon_'+ str(out_sim_pc.lon.values) + '_lat_' + str(out_sim_pc.lat.values) + '.pkl'),"wb"))

            del kp, out_sim_pc

    print('Elapsed time= ' + str(time.time()-start)) 
    

def Reconstruct_KMA_clusters_hindcast(save_path, out_sim, sp_kma, point=[], reconst_hindcast=1):
    '''
    This function generates for each point of the grid, the spectra associated to all of the 2000 KMA clusters and its representative parameters Hs, Tp, Tm, Dir and dspr.
    Also, if reconst_hindcast is activated, reconstructs the hindcast timeseries
    
    save_path          :        Main path for site results folder. Kp: kp_grid folder, efth: efth_kma folder,  stats: reconst_hindcast folder
    out_sim            :        SWAN simulations needed for coordinate locations
    sp_kma             :        dataset with information of the mean energy per cluster
    point              :        If defined, results will only be obtained for the closest point of the grid. If not, all grid points will processed.
    reconst_hindcast   :        If activated, it will reconstruct the hindcast time series 
    '''
    
    if point:
        pos_lon=[np.argmin(np.abs(out_sim.lon.values-point[0]))]
        pos_lat=[np.argmin(np.abs(out_sim.lat.values-point[1]))]
    else:
        pos_lon=range(len(out_sim.lon))
        pos_lat=range(len(out_sim.lat))

    cc=0

    print('Preprocessing and saving efth and stats for all the grid points')
    print('Total: ' + str(len(out_sim.lon) * len(out_sim.lat)))

    a=sp_kma.efth_cluster.expand_dims({'case': np.arange(len(out_sim.case))})
    start=time.time()

    for ilon in pos_lon:
        for ilat in pos_lat:

            sys.stdout.write('\rCount %d' %cc)
            sys.stdout.flush()
            cc+=1
            
            #Reconstruct KMA clusters
            
            if os.path.isfile(os.path.join(save_path, 'efth_kma', 'efth_lon_'+ str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.nc'))== False:
                
                kp= pk.load(open(os.path.join(save_path,  'kp_grid', 'kp_lon_'+ str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.pkl'),'rb'))
                b = np.repeat(kp[:, :, :, np.newaxis], len(sp_kma.clusters), axis=3)
                efth=np.sum(a*b, axis=0)
                efth.to_netcdf(os.path.join(save_path, 'efth_kma', 'efth_lon_'+ str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.nc'))
                #efth['efth_cluster']=(('freq','dir','clusters'),efth.values)
            
            #else:
            efth=xr.open_dataset(os.path.join(save_path, 'efth_kma', 'efth_lon_'+ str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.nc'))                
            
            if os.path.isfile(os.path.join(save_path, 'efth_kma', 'stats_lon_'+ str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.nc')) == False:

                stats = xr.Dataset(
                    {
                        'Hs': (['cluster'], efth.efth_cluster.spec.hs()),
                        'Tp': (['cluster'], efth.efth_cluster.spec.tp()),
                        'Tm': (['cluster'], efth.efth_cluster.spec.tm02()),
                        'Dir_m': (['cluster'], efth.efth_cluster.spec.dm()),
                        'Dir_p': (['cluster'], efth.efth_cluster.spec.dp()),
                        'dspr': (['cluster'], efth.efth_cluster.spec.dspr()),

                    },
                    coords={
                        'cluster': sp_kma.clusters.values,
                    },
                )
                stats.to_netcdf(os.path.join(save_path, 'efth_kma', 'stats_lon_'+ str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.nc'))
                
            else:
                stats=xr.open_dataset(os.path.join(save_path, 'efth_kma', 'stats_lon_'+ str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.nc'))
                          
            #Reconstruct hindcast
            
            if reconst_hindcast:
                
                if os.path.isfile(os.path.join(save_path, 'reconst_hindcast', 'reconst_hind_lon_'+ str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.nc'))==False:

                    output = np.full([len(sp_kma.time), 6], np.nan)
                    for b in range(len(np.unique(sp_kma.bmus))):
                        s=np.where(sp_kma.bmus==b)
                        f=np.array([stats.isel(cluster=b).Hs.values, stats.isel(cluster=b).Tp.values, stats.isel(cluster=b).Tm.values, 
                                    stats.isel(cluster=b).Dir_m.values, stats.isel(cluster=b).Dir_p.values, stats.isel(cluster=b).dspr.values])
                        output[s,:]=np.tile(f,len(s)).reshape(len(s),-1)

                    data = xr.Dataset(
                        {
                            'hs': (['time'], output[:,0]),
                            'tp': (['time'], output[:,1]),
                            'tm': (['time'], output[:,2]),
                            'dm': (['time'], output[:,3]),
                            'dp': (['time'], output[:,4]),
                            'dspr': (['time'], output[:,5]),
                        },
                        coords={
                            'time': sp_kma.time.values,
                        },
                    )     

                    data.to_netcdf(os.path.join(save_path, 'reconst_hindcast', 'reconst_hind_lon_'+ str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.nc'))

    print('Elapsed time= ' + str(time.time()-start)) 
    
    
def SuperPoint_Superposition_forecast(path, date, site, stations,sectors,deg_sup,st_wind):
    
    st=xr.open_dataset(os.path.join(path, 'Data', 'nearshore_forecast',  'spec_forecast_' + date + '_' + site + '_station' + str(stations[0]) + '.nc' ))
        
    efth_all=np.full([len(st.time),len(st.freq),len(st.dir),len(stations)],0.0)
    cont=np.full([len(st.dir)],0)
    del st
    
    for s in range(len(stations)):
        
        print('Station: ' + str(stations[s]))
        
        st=xr.open_dataset(os.path.join(path, 'Data', 'nearshore_forecast',  'spec_forecast_' + date + '_' + site + '_station' + str(stations[s]) + '.nc' ))
        
        if (sectors[s][1]-sectors[s][0])<0:
            d=np.where((st.dir.values>sectors[s][0]-deg_sup) | (st.dir.values<=sectors[s][1]+deg_sup))[0]
        else:
            d=np.where((st.dir.values>sectors[s][0]-deg_sup) & (st.dir.values<=sectors[s][1]+deg_sup))[0]
        
        cont[d] += 1
        efth_all[:,:,d,s]=st.efth[d,:,:].T     
        
        if stations[s]==st_wind:
            wsp=st.Wspeed.values
            wdir=st.Wdir.values
            depth=np.full([len(st.time.values)],st.depth)

    efth_all=(np.sum(efth_all,axis=3)/cont)

    super_point = xr.Dataset({'efth': (['time','freq','dir'],efth_all),
                               'Wspeed': (['time'],wsp),
                               'Wdir': (['time'],wdir),
                               'Depth': (['time'],depth),
                               
                               }, 
                              coords={'time': st.time.values, 
                                      'dir': st.dir.values, 
                                      'freq':st.freq.values})
    super_point['time']=super_point['time'].dt.round('H').values

    return super_point

def Reconstruct_forecast_partitions(path, site, date, sp, out_sim, do=None):

    ini=time.time()

    if do:
        print('Total: ' + str(len(out_sim.lon) * len(out_sim.lat)))

        HS=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)

        HS_SEA=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)
        HS_SW1=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)
        HS_SW2=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)
        HS_SW3=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)

        TP=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)

        TP_SEA=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)
        TP_SW1=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)
        TP_SW2=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)
        TP_SW3=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)   

        TM=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)

        TM_SEA=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)
        TM_SW1=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)
        TM_SW2=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)
        TM_SW3=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)  

        DIR=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)

        DIR_SEA=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)
        DIR_SW1=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)
        DIR_SW2=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)
        DIR_SW3=np.full([len(sp.time.values), len(out_sim.lon), len(out_sim.lat)], np.nan)

        a=sp.efth.expand_dims({'case': np.arange(len(out_sim.case))})

        cc=0

        for ilon in range(len(out_sim.lon)):
            for ilat in range(len(out_sim.lat)):

                sys.stdout.write('\rCount %d' %cc)
                sys.stdout.flush()

                cc+=1

                start=time.time()
                out_sim_p=out_sim.isel(lon=ilon, lat=ilat)

                kp= pk.load(open(os.path.join(path, 'Results', site,  'kp_grid', 'kp_lon_'+ str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.pkl'),'rb'))

                b = np.repeat(kp[:, np.newaxis, :, :], len(sp.time), axis=1)
                sp1=sp.copy()
                sp1['efth']=(('time','freq','dir'),np.sum(a*b, axis=0))
                sp1.to_netcdf(os.path.join(path, 'Results', site, 'nearshore_forecast_efth', 'efth_lon_'+ date + '_'  +  str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.nc'))

                parts=sp1.spec.partition(sp1.Wspeed, sp1.Wdir, sp1.Depth ,wscut=0.00000001,max_swells=3, agefac=1.7).spec.stats(['hs','tp','tm02','dm','dp'])

    #             sys.exit()

                parts.to_netcdf(os.path.join(path, 'Results', site, 'nearshore_forecast_efth', 'parts_lon_'+ date + '_'  +  str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.nc'))

                HS[:,ilon,ilat]=sp1.efth.spec.hs()
                HS_SEA[:,ilon,ilat]=np.where(parts.isel(part=0).hs.values==0, np.nan, parts.isel(part=0).hs.values)
                HS_SW1[:,ilon,ilat]=np.where(parts.isel(part=1).hs.values==0, np.nan, parts.isel(part=1).hs.values)
                HS_SW2[:,ilon,ilat]=np.where(parts.isel(part=2).hs.values==0, np.nan, parts.isel(part=2).hs.values)
                HS_SW3[:,ilon,ilat]=np.where(parts.isel(part=3).hs.values==0, np.nan, parts.isel(part=3).hs.values)


                TP[:,ilon,ilat]=sp1.efth.spec.tp()
                TP_SEA[:,ilon,ilat]=parts.isel(part=0).tp.values
                TP_SW1[:,ilon,ilat]=parts.isel(part=1).tp.values
                TP_SW2[:,ilon,ilat]=parts.isel(part=2).tp.values
                TP_SW3[:,ilon,ilat]=parts.isel(part=3).tp.values

                TM[:,ilon,ilat]=sp1.efth.spec.tm02()
                TM_SEA[:,ilon,ilat]=parts.isel(part=0).tm02.values
                TM_SW1[:,ilon,ilat]=parts.isel(part=1).tm02.values
                TM_SW2[:,ilon,ilat]=parts.isel(part=2).tm02.values
                TM_SW3[:,ilon,ilat]=parts.isel(part=3).tm02.values

                DIR[:,ilon,ilat]=sp1.efth.spec.dm()
                DIR_SEA[:,ilon,ilat]=np.where(parts.isel(part=0).hs.values==0, np.nan, parts.isel(part=0).dm.values)
                DIR_SW1[:,ilon,ilat]=np.where(parts.isel(part=1).hs.values==0, np.nan, parts.isel(part=1).dm.values)
                DIR_SW2[:,ilon,ilat]=np.where(parts.isel(part=2).hs.values==0, np.nan, parts.isel(part=2).dm.values)
                DIR_SW3[:,ilon,ilat]=np.where(parts.isel(part=3).hs.values==0, np.nan, parts.isel(part=3).dm.values)

                del kp, out_sim_p, sp1 #parts, sp1

        print('Elapsed time= ' + str(time.time()-start)) 

        FORECAST = xr.Dataset(
            {
                'hs': (['time','lon','lat'], HS),
                'hs0': (['time','lon','lat'], HS_SEA),
                'hs1': (['time','lon','lat'], HS_SW1),
                'hs2': (['time','lon','lat'], HS_SW2),
                'hs3': (['time','lon','lat'], HS_SW3),

                'tp': (['time','lon','lat'], TP),
                'tp0': (['time','lon','lat'], TP_SEA),
                'tp1': (['time','lon','lat'], TP_SW1),
                'tp2': (['time','lon','lat'], TP_SW2),
                'tp3': (['time','lon','lat'], TP_SW3),

                'tm': (['time','lon','lat'], TM),
                'tm0': (['time','lon','lat'], TM_SEA),
                'tm1': (['time','lon','lat'], TM_SW1),
                'tm2': (['time','lon','lat'], TM_SW2),
                'tm3': (['time','lon','lat'], TM_SW3),

                'dir': (['time','lon','lat'], DIR),
                'dir0': (['time','lon','lat'], DIR_SEA),
                'dir1': (['time','lon','lat'], DIR_SW1),
                'dir2': (['time','lon','lat'], DIR_SW2),
                'dir3': (['time','lon','lat'], DIR_SW3),
            },
            coords={
                'time': sp.time.values,
                'lon' : out_sim.lon,
                'lat' : out_sim.lat
            },
        )
        FORECAST.to_netcdf(os.path.join(path, 'Results', site, 'nearshore_forecast', 'Forecast_grid_bulk_parts' + date + '_' + site + '.nc' ))

    else:
        FORECAST=xr.open_dataset(os.path.join(path, 'Results', site, 'nearshore_forecast', 'Forecast_grid_bulk_parts' + date + '_' + site + '.nc' ))

    return FORECAST


def Get_hindcast_time(path, site, out_sim, time_sel):

    if os.path.isfile(os.path.join(path, 'Results', site, 'reconst_hindcast','Hindcast_time_'+ time_sel +'.nc'))==False:

        pos_lon=range(len(out_sim.lon))
        pos_lat=range(len(out_sim.lat))
        cc=0

        hs_metrics=np.full([len(pos_lon), len(pos_lat)], np.nan)
        tp_metrics=np.full([len(pos_lon), len(pos_lat)], np.nan)
        tm_metrics=np.full([len(pos_lon), len(pos_lat)], np.nan)
        dir_metrics=np.full([len(pos_lon), len(pos_lat)], np.nan)

        DATA_S=[]

        for ilon in pos_lon:
                for ilat in pos_lat:

                    file_load=os.path.join(path, 'Results', site, 'reconst_hindcast', 'reconst_hind_lon_'+ str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.nc')

                    if os.path.isfile(file_load)==True:

                        sys.stdout.write('\rCount %d' %cc)
                        sys.stdout.flush()
                        cc+=1

                        data=xr.open_dataset(file_load).sel(time=time_sel)

                        hs_metrics[ilon, ilat] = data.hs.values
                        tp_metrics[ilon, ilat] = data.tp.values
                        tm_metrics[ilon, ilat] = data.tm.values
                        dir_metrics[ilon, ilat] = data.dm.values

        time_stats = xr.Dataset({'Hs': (['lon','lat'],hs_metrics),
                                 'Tp': (['lon','lat'],tp_metrics),
                                 'Tm': (['lon','lat'],tm_metrics),
                                 'Dir': (['lon','lat'],dir_metrics),
                                 'Fe': (['lon','lat'],(hs_metrics**2)*tp_metrics),

                   }, 
                  coords={'lon': out_sim.lon.values,
                          'lat': out_sim.lat.values, 
                         })
        time_stats.to_netcdf(os.path.join(path, 'Results', site, 'reconst_hindcast','Hindcast_time_'+ time_sel +'.nc'))
    else:
        
        time_stats=xr.open_dataset(os.path.join(path, 'Results', site, 'reconst_hindcast','Hindcast_time_'+ time_sel +'.nc'))
        
        
    return time_stats
    
    
def obtain_wave_metrics_grid(path, site, out_sim, perc=[25,50,75,90,95,99]):
    
    '''
    Obtain general metrics for plotting
    Generates mean direction for all the points in the grid and hs, tp and tm for the different percentiles defined in perc
    '''
    
    pos_lon=range(len(out_sim.lon))
    pos_lat=range(len(out_sim.lat))
    cc=0

    hs_metrics=np.full([len(pos_lon), len(pos_lat), len(perc)], np.nan)
    tp_metrics=np.full([len(pos_lon), len(pos_lat), len(perc)], np.nan)
    tm_metrics=np.full([len(pos_lon), len(pos_lat), len(perc)], np.nan)
    dir_mean=np.full([len(pos_lon), len(pos_lat)], np.nan)
    dir_peak=np.full([len(pos_lon), len(pos_lat)], np.nan)

    DATA_S=[]

    for ilon in pos_lon:
            for ilat in pos_lat:

                file_load=os.path.join(path, 'Results', site, 'reconst_hindcast', 'reconst_hind_lon_'+ str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.nc')

                if os.path.isfile(file_load)==True:

                    sys.stdout.write('\rCount %d' %cc)
                    sys.stdout.flush()

                    data=xr.open_dataset(file_load)

                    hs_metrics[ilon, ilat, :]=np.nanpercentile(data.hs.values, perc)
                    tp_metrics[ilon, ilat, :]=np.nanpercentile(data.tp.values, perc)
                    tm_metrics[ilon, ilat, :]=np.nanpercentile(data.tm.values, perc)

                    dir_mean[ilon,ilat]=circmean(data.dm.values, high = 360, low = 0)
                    dir_peak[ilon,ilat]=circmean(data.dp.values, high = 360, low = 0)
                    cc+=1

                    data_s = data.resample(time='1M').mean().assign_coords(lon=out_sim.lon.values[ilon], lat=out_sim.lat.values[ilat]).expand_dims(['lon','lat'])

                    if not DATA_S:
                        DATA_S=data_s
                    else:
                        DATA_S=DATA_S.merge(data_s)
    
    hindcast_stats = xr.Dataset({'Hs': (['lon','lat','perc'],hs_metrics),
                             'Tp': (['lon','lat','perc'],tp_metrics),
                             'Tm': (['lon','lat','perc'],tm_metrics),
                             'Dir_m': (['lon','lat'],dir_mean),
                             'Dir_p': (['lon','lat'],dir_peak),
                             'Fe': (['lon','lat','perc'],(hs_metrics**2)*tp_metrics),

               }, 
              coords={'lon': out_sim.lon.values,
                      'lat': out_sim.lat.values, 
                      'perc': perc,})
    
    return hindcast_stats, DATA_S

    
#############################################################
########################   PLOTTING  ########################
#############################################################


import random
def Plot_obstructions_CAWCR(obstr_cawcr, fname, area=[], min_z=-10000, max_z=1000, figsize=[25,12]):
    
    '''
    obstr_cawcr  : xarray with obstructions from CAWCR
    area         : if defined, map box
    
     '''
    
    coast = ShapelyFeature(Reader(fname).geometries(),
                                    ccrs.PlateCarree(), edgecolor='black')
    
    fig = plt.figure(random.randint(1,50), figsize=figsize)
    gs3=gridspec.GridSpec(1,1)
    ax=fig.add_subplot(gs3[0],projection = ccrs.PlateCarree(central_longitude=180))

    ax.add_feature(coast,facecolor="None", alpha=0.8, edgecolor='black', linewidth=1.5, zorder=5)
    
    if area:
        pos_lon=np.where((obstr_cawcr.lon.values>area[0]) & (obstr_cawcr.lon.values<area[1]))[0]
        pos_lat=np.where((obstr_cawcr.lat.values>area[2]) & (obstr_cawcr.lat.values<area[3]))[0]
        obstr_cawcr=obstr_cawcr.isel(lon=pos_lon, lat=pos_lat)
        ax.set_extent([area[0], area[1], area[2], area[3]], crs=ccrs.PlateCarree())

    c1=colormap(obstr_cawcr.depth.values/1000, max_z, min_z)
    im=ax.pcolormesh(obstr_cawcr.lon, obstr_cawcr.lat, obstr_cawcr.depth/1000, cmap = c1, vmax = max_z, vmin = min_z,transform=ccrs.PlateCarree(),zorder=1, alpha=0.3)

    obstr_cawcr.obstructions_y.values[np.where(obstr_cawcr.obstructions_y.values==0)]=np.nan
    obstr_cawcr.obstructions_x.values[np.where(obstr_cawcr.obstructions_x.values==0)]=np.nan
    
    xx,yy=np.meshgrid(obstr_cawcr.lon,obstr_cawcr.lat)
    sc_y=ax.pcolormesh(obstr_cawcr.lon, obstr_cawcr.lat, obstr_cawcr.obstructions_y, cmap='Purples', vmin=0, vmax=100, transform=ccrs.PlateCarree(), zorder=2, alpha=0.6)
    sc_x=ax.pcolormesh(obstr_cawcr.lon, obstr_cawcr.lat, obstr_cawcr.obstructions_x, cmap='Reds', vmin=0, vmax=100, transform=ccrs.PlateCarree(), zorder=2, alpha=0.6)
    
    if not area:
        ax.contour(obstr_cawcr.lon, obstr_cawcr.lat, obstr_cawcr.mask, 3, cmap='viridis', zorder=3,transform=ccrs.PlateCarree())
        
    if area:
        ax.set_extent([area[0], area[1], area[2], area[3]], crs=ccrs.PlateCarree())

    # plt.colorbar(im).set_label('Elevation (m)')
    plt.colorbar(sc_y).set_label('Obstruction Y (%)', fontsize=16)
    plt.colorbar(sc_x).set_label('Obstruction X (%)', fontsize=16)


def colormap(elev, topat, topag): #Elevation, maximum topo, minimum bati
    colors2 = 'YlGnBu_r'
    colors1 = cmocean.cm.turbid
    bottom = plt.get_cmap(colors2, -topag)
    top = plt.get_cmap(colors1, topat)
    newcolors = np.vstack((bottom(np.linspace(0, 0.8, -topag)),top(np.linspace(0.1, 1, topat))))
    return(ListedColormap(newcolors))

def plot_bati(bati,fname=[], min_z=-200, max_z=15, figsize=[25,15]):
    
    fig = plt.figure(figsize=figsize)
    gs3=gridspec.GridSpec(1,1)
    ax=fig.add_subplot(gs3[0],projection = ccrs.PlateCarree(central_longitude=180))
    
    if fname:
        shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                    ccrs.PlateCarree(), edgecolor='black')

        ax.add_feature(shape_feature,facecolor="None", alpha=1, edgecolor='black', linewidth=1.5, zorder=2)
        
    c1=colormap(-bati.z, max_z, min_z)
    im=ax.pcolormesh(bati.x.values, bati.y.values,-bati.z.values, cmap = c1, vmax = max_z, vmin = min_z,transform=ccrs.PlateCarree(),zorder=1 )
    plt.colorbar(im)    .set_label('Elevation (m)')
    
def plot_bati_3D(bati, min_z=-200, max_z=15):
    
    fig=plt.figure(figsize=[17,14])
    c1=colormap(-bati.z, max_z, min_z)
    ax = fig.add_subplot(111, projection='3d') 
    xx,yy=np.meshgrid(bati.x.values, bati.y.values)
    ax.plot_surface(xx,yy, -bati.z,cmap = c1, vmax = max_z, vmin = min_z)


def plot_spectrum(ax,x,y,z, vmax=0.3, ylim=0.49, cmap='magma_r', vmin=0, remove_axis=0):
    
    x1=np.append(x,x[0])
    y1=np.append(0,y)
    z1=np.column_stack((z[:,:],z[:,-1]))
    
    if cmap=='RdBu_r':
        norm = mcolors.TwoSlopeNorm(0,vmin, vmax)
        p1=ax.pcolormesh(x1,y1,z1, vmin=vmin,vmax=vmax, cmap=plt.cm.RdBu_r, norm=norm)
    else:
        p1=ax.pcolormesh(x1,y1,z1, vmin=vmin, vmax=vmax)   
    
    p1.set_cmap(cmap)
    ax.set_theta_zero_location('N', offset=0)
    ax.set_theta_direction(-1)
    ax.set_ylim(0,ylim)
    
    if remove_axis:
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        ax.tick_params(axis='y', colors='plum',labelsize=14,grid_linestyle=':',grid_alpha=0.75,grid_color='plum')
        ax.tick_params(axis='x', colors='purple',labelsize=14,pad=5,grid_linestyle=':',grid_alpha=0.75,grid_color='plum')
    
    ax.grid(color='plum', linestyle='--', linewidth=0.7,alpha=0.2)
    return p1

def Plot_empty_grid(sp, pos=[25,6]):
    
    '''
    Function to exemplify the binwaves process of selecting individual grid cells
    
    sp   :    spectra
    pos  :    position to highlight grid
    '''
    
    fig = plt.figure(figsize=[8,8])
    gs3=gridspec.GridSpec(1,1,hspace=0.01, wspace=0.01) 
    ax=fig.add_subplot(gs3[0],projection='polar')

    x1=np.append(np.deg2rad(sp.dir.values-7.5),np.deg2rad(sp.dir.values-7.5)[0])
    y1=np.append(0,sp.freq.values)
    z=np.full(np.shape(sp.efth.values[0,:,:]),0.0)
    z[pos[0],pos[1]]=1
    z1=np.column_stack((z[:,:],z[:,-1]))

    p1=ax.pcolormesh(x1,y1,z1, vmin=0, vmax=1, edgecolor='powderblue', linewidth=0.05, cmap='Reds')   
    ax.set_theta_zero_location('N', offset=0)
    ax.set_theta_direction(-1)
    ax.tick_params(axis='both', colors='royalblue',labelsize=17)
    
def Plot_cases_grid(sp, cases_id):
    
    '''
    Function to plot the cases with different colors
    
    sp   :    spectra
    cases_id  :    matrix with case identifyer for each freq and dir
    '''
    
    fig = plt.figure(figsize=[10,10])
    gs3=gridspec.GridSpec(1,1,hspace=0.01, wspace=0.01) 
    ax=fig.add_subplot(gs3[0],projection='polar')

    x1=np.append(np.deg2rad(sp.dir.values-7.5),np.deg2rad(sp.dir.values-7.5)[0])
    y1=np.append(0,sp.freq.values)
    z=cases_id
    
    z1=np.column_stack((z[:,:],z[:,-1]))
        
    cmn = np.vstack((cm.get_cmap('plasma', 124)(np.linspace(0, 0.9, 70)),
                     cm.get_cmap('magma_r', 124)(np.linspace(0.1, 0.4, 80)),
                     cm.get_cmap('rainbow_r', 124)(np.linspace(0.1, 0.8, 80)),
                     cm.get_cmap('Blues_r', 124)(np.linspace(0.4, 0.8, 40)),
                     cm.get_cmap('cubehelix_r', 124)(np.linspace(0.1, 0.8, 80)),
                     ))
    
    cmn = ListedColormap(cmn, name='cmn')

    p1=ax.pcolormesh(x1,y1,z1, vmin=0, vmax=np.nanmax(cases_id), edgecolor='grey', linewidth=0.01, cmap=cmn, shading='flat')   
    ax.set_theta_zero_location('N', offset=0)
    ax.set_theta_direction(-1)
    ax.tick_params(axis='both', colors='black',labelsize=14, pad=10)
    plt.colorbar(p1, pad=0.1, shrink=0.7).set_label('Case ID', fontsize=16)
    
def plot_swan_case(out_sim,subset, case,  lon_s=[], lat_s=[], coast=[], var_list=None):
    
    fig = plt.figure(figsize=[27,7])
    if not var_list:
        var_list=['Hsig','TPsmoo','Dir','Dspr']
        
    gs3=gridspec.GridSpec(1,len(var_list), hspace=0.01, wspace=0.01)

    for a in range(len(var_list)):

        if var_list[a] =='Dir': 
            cmap='twilight_shifted'
            vmin, vmax=0, 360
        if var_list[a] =='Hsig': 
            cmap='RdBu_r'
            vmin, vmax=0, 2
        if var_list[a] =='TPsmoo': 
            cmap='magma_r'
            vmin, vmax=0, 30
        if var_list[a] =='Dspr': 
            cmap='rainbow'
            vmin, vmax=0, 80

        ax=fig.add_subplot(gs3[a],projection = ccrs.PlateCarree(central_longitude=180))
        if coast:
            ax.add_feature(coast,facecolor="None", alpha=0.5, edgecolor='navy', linewidth=1.5, zorder=2)
        ax.gridlines()
        XX, YY = np.meshgrid(out_sim.lon.values, out_sim.lat.values)
        pm = ax.pcolormesh(
                XX,YY, out_sim[var_list[a]].sel(case=case).values.T,
                cmap=cmap,
                vmin=vmin, vmax=vmax,
                shading='auto',
                transform=ccrs.PlateCarree()
        )
        #ax.set_xlim([np.nanmin(out_sim.lon.values), np.nanmax(out_sim.lon.values)])
        #ax.set_ylim([np.nanmin(out_sim.lat.values), np.nanmax(out_sim.lat.values)])
        plt.colorbar(pm, orientation="horizontal", pad=0.1).set_label(var_list[a])
        if lat_s:
            ax.plot(lon_s, lat_s, 'k*', transform=ccrs.PlateCarree(), markersize=15)
        
        ax.set_title(var_list[a], fontsize=15)
    plt.suptitle('Tp= ' + str(subset.sel(case=case).per.values) + '   Dir= ' + str(subset.sel(case=case).dir.values), fontsize=16, y=0.98)
    
def plot_spec_prop(input_spec, output_spec,sp, c, ylim=0.5, cmap='gnuplot2_r'):
    
    fig = plt.figure(figsize=[25,7])
    var_list=['Hsig','TPsmoo','Dir','Dspr']
    gs3=gridspec.GridSpec(1,2, hspace=0.01, wspace=0.01)
    
    ax=fig.add_subplot(gs3[0],projection = 'polar')
    ax.set_title('Input', fontsize=16)
    y=sp.freq.values
    x= np.deg2rad(sp.dir.values)
    z=input_spec.sel(case=c).hs.values
    p1=plot_spectrum(ax,x ,y, z, vmax=1, ylim=ylim, cmap=cmap)
    plt.colorbar(p1).set_label('Hs')
    
    ax=fig.add_subplot(gs3[1],projection = 'polar')
    ax.set_title('Propagated to Point', fontsize=16)
    y=sp.freq.values
    x= np.deg2rad(sp.dir.values)
    z=output_spec.sel(case=c).kp.values
    p1=plot_spectrum(ax,x ,y, z, vmax=1, ylim=ylim, cmap=cmap)
    plt.colorbar(p1).set_label('Kp')
    
    
def plot_reconstructed_spec(sp_in, s_p ,sp, vmax=0.5, ylim=0.5, cmap='gnuplot2_r'):
    
    hs_in=sp_in.efth.spec.hs().values
    hs_p=sp_in.efth_p.spec.hs().values
    
    fig = plt.figure(figsize=[25,7])
    var_list=['Hsig','TPsmoo','Dir','Dspr']
    gs3=gridspec.GridSpec(1,2, hspace=0.01, wspace=0.01)
    
    ax=fig.add_subplot(gs3[0],projection = 'polar')
    ax.set_title('Input Hs= ' + str(np.round(hs_in,3)), fontsize=16)
    y=sp.freq.values
    x= np.deg2rad(sp.dir.values)
    z=np.sqrt(sp_in.efth.values)
    p1=plot_spectrum(ax,x ,y, z, vmax=vmax, ylim=ylim, cmap=cmap)
    plt.colorbar(p1).set_label('Sqrt(Efth)')
    
    ax=fig.add_subplot(gs3[1],projection = 'polar')
    ax.set_title('Propagated Hs= ' + str(np.round(hs_p,3)), fontsize=16)
    y=sp.freq.values
    x= np.deg2rad(sp.dir.values)
    z=np.sqrt(s_p)
    p1=plot_spectrum(ax,x ,y, z, vmax=vmax, ylim=ylim, cmap=cmap)
    plt.colorbar(p1).set_label('Sqrt(Efth)')

    
def plot_EOFS_spec_energy(sp_mod_daily, gr=5):
    
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
    
    
def Plot_kmeans_clusters(sp, annomaly=[], vmax=1.2, num_clusters=[], figsize=[16,12], cmap='CMRmap_r', ylim=0.5):
    
    if not num_clusters:
        num_clusters=len(np.unique(sp.bmus))

    mean_z=np.sqrt(np.mean(sp.efth.values,axis=0))
    
    fig = plt.figure(figsize=figsize)

    h_clust=np.full([len(sp.freq), len(sp.dir), len(sp.kma_order)], np.full)

    gs3=gridspec.GridSpec(np.int(np.sqrt(num_clusters)),np.int(np.sqrt(num_clusters)),hspace=0.01, wspace=0.01)
    hs_m=[]
    for b in range(num_clusters):
        s=np.where(sp.bmus==b)[0]
        prob=len(s)/len(sp.bmus)
        ax2=fig.add_subplot(gs3[b],projection='polar')
        if annomaly:
            z=np.sqrt(np.mean(sp.efth.isel(time=s).values,axis=0)) - mean_z
            h_clust[:,:,b]=z
            y=sp.freq.values
            x= np.deg2rad(sp.dir.values)
            p2=plot_spectrum(ax2,x,y,z,vmin=-vmax, vmax=vmax,ylim=np.nanmax(sp.freq), cmap='RdBu_r', remove_axis=1)  
        else:
            z=np.sqrt(np.mean(sp.efth.isel(time=s).values,axis=0))
            y=sp.freq.values
            x= np.deg2rad(sp.dir.values)
            p2=plot_spectrum(ax2,x,y, z, vmin=0, vmax=vmax, ylim=np.nanmax(sp.freq), remove_axis=1, cmap=cmap) 
            h_clust[:,:,b]=z
        ax2.set_ylim([0, ylim])
        hs_m.append(np.sqrt(np.sum(z**2)))
        
        
    gs3.tight_layout(fig, rect=[[], [], 0.83, []])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p2,cax=ax0)
    if annomaly:
        ax0.set_ylabel('sqrt(Efth) - ANNOMALIES', fontsize=14)
    else:
        ax0.set_ylabel('sqrt(Efth)', fontsize=14)
    gs4.tight_layout(fig, rect=[0.83, 0.1, 0.9, 0.9])
    
    return h_clust, hs_m

def plot_specs_inside_cluster(sp, cluster, gr=6, vmax=0.7, ylim=0.3):
    
    fig = plt.figure(figsize=[14,12])

    s=np.where(sp.bmus==cluster)[0]

    gs3=gridspec.GridSpec(gr,gr,hspace=0.01, wspace=0.01)
    hs_m=[]

    if len(s)<gr**2:
        lim=len(s)
    else:
        lim=gr**2

    for b in range(np.int(lim)):

        ax2=fig.add_subplot(gs3[b],projection='polar')
        z=np.sqrt(sp.efth.isel(time=s[b]).values)
        y=sp.freq.values
        x= np.deg2rad(sp.dir.values)
        hs_m.append(np.sqrt(np.sum(z**2)))
        cmap='magma_r'
        p2=plot_spectrum(ax2,x,y, z, vmin=0, vmax=vmax, ylim=np.nanmax(sp.freq), remove_axis=1, cmap=cmap)
        ax2.set_ylim([0, ylim])

        gs3.tight_layout(fig, rect=[[], [], 0.9, []])

    gs4=gridspec.GridSpec(1,1)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p2,cax=ax0)
    ax0.set_ylabel('Hs (m)', fontsize=14)
    gs4.tight_layout(fig, rect=[0.9, 0.1, 0.99, 0.9])
    
    
def example_plot_spectra(save_path, out_sim, point, sp_kma, cluster):
    
    '''
    Plot of example spectra from a KMA cluster in a specific point
    '''
    
    ilon=np.argmin(np.abs(out_sim.lon.values-point[0]))
    ilat=np.argmin(np.abs(out_sim.lat.values-point[1]))

    efth=xr.open_dataset(os.path.join(save_path, 'efth_kma', 'efth_lon_'+ str(out_sim.lon.values[ilon]) + '_lat_' + str(out_sim.lat.values[ilat]) + '.nc'))

    fig = plt.figure(figsize=[8,8])
    gs3=gridspec.GridSpec(1,1,hspace=0.01, wspace=0.01)
    y=sp_kma.freq.values
    x= np.deg2rad(sp_kma.dir.values)
    z=np.sqrt(efth.efth_cluster.values[:,:,cluster])
    ax2=fig.add_subplot(gs3[0],projection='polar')
    vmax=0.
    p2=plot_spectrum(ax2,x,y, z, vmin=0, vmax=0.8, ylim=np.nanmax(sp_kma.freq), remove_axis=1, cmap='CMRmap_r') 
    
    
    
def Plot_comparison_buoy(data, buoy_data, ss=3):
    '''
    This function plots the comparison of the reconstruction against buoy data
    '''
    
    fig = plt.figure(1, figsize=[20,14])
    #Timeseries
    gs1=gridspec.GridSpec(3,1)
    ax1=fig.add_subplot(gs1[0])
    ax2=fig.add_subplot(gs1[1],sharex=ax1)
    ax3=fig.add_subplot(gs1[2],sharex=ax1)
    color=['navy', 'crimson' ,'darkmagenta','springgreen','purple','lightseagreen','indianred','orange','orchid' ]

    ax1.plot(data.time,data.hs,'.-',markersize=ss, color='darkmagenta', label='BinWaves')
    ax1.plot(buoy_data.time,buoy_data.hs,'.-',markersize=ss, color='navy', label='buoy')
    ax1.set_ylim([0, 1.1*np.nanmax(buoy_data.hs)])
    ax1.set_ylabel('Hs (m)')
    ax1.legend()

    ax2.plot(data.time,data.tp,'.-', markersize=ss, color='plum', label='Tp - BinWaves')
    
    
    try:
        np.nanmax(buoy_data.tp)
        ax2.plot(buoy_data.time,buoy_data.tp,'.-', markersize=ss, color='royalblue', label='Tp - Buoy')
    except:
        pass
    
    try:
        np.nanmax(buoy_data.tm02)
        ax2.plot(buoy_data.time,buoy_data.tm02,'.-', markersize=ss, color='navy', label='T02 - Buoy')
        ax2.plot(data.time,data.tm,'.-', markersize=ss, color='magenta', label='Tm - BinWaves')
    except:
        pass
    
#     try:
#         np.nanmax(buoy_data.tm10)
#         ax2.plot(buoy_data.time,buoy_data.tm10,'.-', markersize=ss, color='skyblue', label='T10 - Buoy')
#     except:
#         pass
    
    ax2.legend(ncol=5)
    ax2.set_ylabel('T (s)')

    ax3.plot(data.time,data.dm,'.',markersize=ss, color='crimson')
    ax3.set_ylabel('Dpm (º)')
    ax3.set_xlim([buoy_data.time.values[0],buoy_data.time.values[-1]])
    
    gs1.tight_layout(fig, rect=[0.01, 0.1, 0.7, 0.9])
    
    #Scatter
    gs2=gridspec.GridSpec(1,1)
    ax1=fig.add_subplot(gs2[0])

    ax1.plot([0,np.nanmax(buoy_data.hs)], [0,np.nanmax(buoy_data.hs)], color='darkred')
    c,ia,ib = np.intersect1d(data.time, buoy_data.time, return_indices=True)
    xy = np.vstack([buoy_data.hs[ib],data.hs[ia]])
    z = gaussian_kde(xy)(xy)

    ax1.scatter(buoy_data.hs[ib],data.hs[ia], c=z, s=50)
    ax1.set_ylabel('Hs - BinWaves KMA', fontsize=16)
    ax1.set_xlabel('Hs - Buoy', fontsize=16)
    ax1.set_xlim([0,np.nanmax(buoy_data.hs)])
    ax1.set_ylim([0,np.nanmax(buoy_data.hs)])
    ax1.set_aspect('equal', 'box')
    
    gs2.tight_layout(fig, rect=[0.75, 0.3, 1, 0.7])
    
def plot_hindcast_series_point(h_p, data_p, var):
    
    '''
    Plot time series of var in the hindcast period
    '''
    
    fig = plt.figure(figsize=[23,7])
    gs3=gridspec.GridSpec(1,1,hspace=0.01, wspace=0.01) 
    ax=fig.add_subplot(gs3[0])

    if (var=='dm') | (var=='dp'):
        ax.plot(h_p.time, h_p[var], '.', color='royalblue', markersize=1, alpha=0.3, label='hourly')
        ax.plot(data_p.time, data_p[var],'.', color='darkmagenta', markersize=10,  label='monthly')
        ax.plot(data_p.resample(time='1Y').mean().time, data_p.resample(time='1Y').mean()[var],'.', markersize=20, color='navy', alpha=1, label='annual')
    else:
        ax.plot(h_p.time, h_p[var], color='royalblue', alpha=0.4, label='hourly')
        ax.plot(data_p.time, data_p[var], color='darkmagenta', alpha=0.7, linewidth=1.4, label='monthly')
        ax.plot(data_p.resample(time='1Y').mean().time, data_p.resample(time='1Y').mean()[var], linewidth=1.5,color='navy', alpha=0.7, label='annual')
    ax.set_xlabel('Time', fontsize=18, fontweight='bold')
    ax.set_ylabel(var, fontsize=18, fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlim([h_p.time.values[0], h_p.time.values[-1]])
    ax.legend(fontsize=16)    
    
def Plot_hindcast_maps(hindcast_stats, coast, var_s, p=[], vmin=[], vmax=[], figsize=[15,15], ax=[], quiv_fact=10, quiv_color='lightslategrey'):
    '''
    Function to plot hindcast maps
    '''
    if p:
        pi=np.argmin(np.abs(hindcast_stats.perc.values-np.array(p)))
    else:
        pi=0
    
    if not vmin:
        vmin=0.9*np.nanmin(hindcast_stats[var_s])
        
    if not vmax:
        vmax=1.1*np.nanmax(hindcast_stats[var_s])
    
    if not ax:
        
        fig = plt.figure(figsize=figsize)
        gs3=gridspec.GridSpec(1,1, hspace=0.01, wspace=0.01)
        ax=fig.add_subplot(gs3[0],projection = ccrs.PlateCarree(central_longitude=180))

    ax.add_feature(coast,facecolor="grey", alpha=1, edgecolor='grey', linewidth=1.5, zorder=20)
    ax.gridlines()

    XX, YY = np.meshgrid(hindcast_stats.lon.values, hindcast_stats.lat.values)

    if (var_s=='Tm') | (var_s=='Tp') | (var_s=='tm') | (var_s=='tp'):
        delta=0.2
        cmn='viridis_r'    

    if (var_s=='Hs') | (var_s=='hs'):
        delta=0.05
        cmn = np.vstack((cm.get_cmap('bwr_r', 124)(np.linspace(0.5, 0.97, 20)),
                         cm.get_cmap('jet', 124)(np.linspace(0.2, 0.84, 80)),
                         ))
        cmn = ListedColormap(cmn, name='cmn')
    if (var_s=='Fe'):
        delta=5
        cmn='plasma_r'
        
    if (var_s=='Dir') | (var_s=='dir') | (var_s=='Dir_m') | (var_s=='Dir_p')  :
        delta=1
        cmn='twilight_shifted'
        
    if not p:
        
        levels=np.arange(0.9*vmin, 1.1*vmax, delta)
        pm = ax.contourf(
                XX,YY, hindcast_stats[var_s].values.T, levels,
                cmap=cmn,
                vmin=vmin, vmax=vmax,
                transform=ccrs.PlateCarree(),
                zorder=2, alpha=0.8
        )
        
    else:

        levels=np.arange(0.9*vmin, 1.1*vmax, delta)
        pm = ax.contourf(
                XX,YY, hindcast_stats[var_s].values[:,:,pi].T, levels,
                cmap=cmn,
                vmin=vmin, vmax=vmax,
                transform=ccrs.PlateCarree(),
                zorder=2, alpha=0.9
        )
        
    try:
        
        flon=np.arange(np.floor(quiv_fact/2).astype('int'),len(hindcast_stats.lon.values), quiv_fact)
        flat=np.arange(np.floor(quiv_fact/2).astype('int'),len(hindcast_stats.lat.values), quiv_fact)
        mask=np.meshgrid(flon,flat)
        
#         ax.streamplot(XX,YY, (np.sin(np.deg2rad(hindcast_stats.Dir_m.values-180)).T),
#                        (np.cos(np.deg2rad(hindcast_stats.Dir_m.values.T-180))), 
#                        zorder=15, linewidth=1.2, color='mediumslateblue',  transform=ccrs.PlateCarree(),
#                        arrowsize=2, density=1.2
#                        )
        hindcast_stats['dm']=(('lon','lat'),np.where(hindcast_stats['dm'].values==270, np.nan, hindcast_stats['dm'].values))
    
        ax.quiver(hindcast_stats.lon.values[flon], 
               hindcast_stats.lat.values[flat], 
               (np.sin(np.deg2rad(hindcast_stats.dm.values[mask]-180))),
               (np.cos(np.deg2rad(hindcast_stats.dm.values[mask]-180))),
               #hindcast_stats.Hs.values[:,:,pi][mask], cmap=cmn, clim=(vmin,vmax),
               color=quiv_color, zorder=15, alpha=0.6,  transform=ccrs.PlateCarree(), width=0.003 , headwidth=5, label='Sea'
              )
        
        
    except:
               
        pass
        
#         ax.streamplot(XX,YY, (np.sin(np.deg2rad(hindcast_stats.dm.values-180)).T),
#                        (np.cos(np.deg2rad(hindcast_stats.dm.values.T-180))), 
#                        zorder=15, linewidth=1.2, color='mediumslateblue',  transform=ccrs.PlateCarree(),
#                        arrowsize=2, density=1.2
#                        )
    try:
        
        flon=np.arange(np.floor(quiv_fact/2).astype('int'),len(hindcast_stats.lon.values), quiv_fact)
        flat=np.arange(np.floor(quiv_fact/2).astype('int'),len(hindcast_stats.lat.values), quiv_fact)
        mask=np.meshgrid(flon,flat)
        
        hindcast_stats['Dir_m']=(('lon','lat'),np.where(hindcast_stats['Dir_m'].values==270, np.nan, hindcast_stats['Dir_m'].values))
                
        ax.quiver(hindcast_stats.lon.values[flon], 
               hindcast_stats.lat.values[flat], 
               (np.sin(np.deg2rad(hindcast_stats.Dir_m.values[mask]-180))),
               (np.cos(np.deg2rad(hindcast_stats.Dir_m.values[mask]-180))),
               #hindcast_stats.Hs.values[:,:,pi][mask], cmap=cmn, clim=(vmin,vmax),
               color=quiv_color, zorder=15, alpha=0.6, transform=ccrs.PlateCarree(), width=0.002, headwidth=6, label='Sea'
              )
    except:
        pass
    
    try:
        
        flon=np.arange(np.floor(quiv_fact/2).astype('int'),len(hindcast_stats.lon.values), quiv_fact)
        flat=np.arange(np.floor(quiv_fact/2).astype('int'),len(hindcast_stats.lat.values), quiv_fact)
        mask=np.meshgrid(flon,flat)
        
#         ax.streamplot(XX,YY, (np.sin(np.deg2rad(hindcast_stats.Dir_m.values-180)).T),
#                        (np.cos(np.deg2rad(hindcast_stats.Dir_m.values.T-180))), 
#                        zorder=15, linewidth=1.2, color='mediumslateblue',  transform=ccrs.PlateCarree(),
#                        arrowsize=2, density=1.2
#                        )
        hindcast_stats['dm']=(('lon','lat'),np.where(hindcast_stats['Dir'].values==270, np.nan, hindcast_stats['dm'].values))
    
        ax.quiver(hindcast_stats.lon.values[flon], 
               hindcast_stats.lat.values[flat], 
               (np.sin(np.deg2rad(hindcast_stats.dm.values[mask]-180))),
               (np.cos(np.deg2rad(hindcast_stats.dm.values[mask]-180))),
               #hindcast_stats.Hs.values[:,:,pi][mask], cmap=cmn, clim=(vmin,vmax),
               color=quiv_color, zorder=15, alpha=0.6,  transform=ccrs.PlateCarree(), width=0.003 , headwidth=5, label='Sea'
              )
    except:
        pass
        


    plt.colorbar(pm, orientation="horizontal", pad=0.1).set_label(var_s, fontsize=16)
    if p:
        ax.set_title(var_s + ' (perc = ' + str(p) + ')', fontsize=20, pad=25)
    else:
        ax.set_title(var_s, fontsize=20, pad=25)
        
    ax.set_extent([np.nanmin(hindcast_stats.lon.values), np.nanmax(hindcast_stats.lon.values), np.nanmin(hindcast_stats.lat.values), np.nanmax(hindcast_stats.lat.values)])
        
    
def Plot_Forecast_Time(FORECAST, time_plot, coast, sp, out_sim, vmin=[], vmax=[], quiv_fact=1, figsize=[29,11]):
    
    if not vmin:
        vmin=np.nanmin(FORECAST.hs.values)
    if not vmax:
        vmax=np.nanmax(FORECAST.hs.values)
    
    fig = plt.figure(figsize=figsize)
    gs3=gridspec.GridSpec(1,6,hspace=0.01, wspace=0.01)

    ax1=fig.add_subplot(gs3[0:2],projection='polar')
    z=sp.isel(time=time_plot).efth.values
    plot_spectrum(ax1,np.deg2rad(sp.dir.values), sp.freq.values,z, cmap='CMRmap_r', ylim=0.3)
    ax1.set_title('Super Point \n ' + str(sp.isel(time=time_plot).time.values.astype('datetime64[h]')), fontsize=16, color='darkmagenta')

    ax2=fig.add_subplot(gs3[2:6],projection = ccrs.PlateCarree(central_longitude=180))
    ax2.add_feature(coast,facecolor='grey', alpha=1, edgecolor='grey', linewidth=1.5, zorder=4)

    cmn = np.vstack((cm.get_cmap('bwr_r', 124)(np.linspace(0.5, 0.96, 20)),
                             cm.get_cmap('jet', 124)(np.linspace(0.2, 0.84, 80)),
                             ))
    cmn = ListedColormap(cmn, name='cmn')

    im=ax2.pcolormesh(FORECAST.lon.values, FORECAST.lat.values, FORECAST.isel(time=time_plot).hs.T, vmin=vmin, vmax=vmax, cmap=cmn, transform=ccrs.PlateCarree(),zorder=1, 
                      alpha=0.4, edgecolor='face',
                      antialiased=True, linewidth=0.15)
    
    cbar=ax2.figure.colorbar(im, ax=ax2)
    cbar.set_label('Hs (m)', fontsize=18)
    cbar.ax.tick_params(labelsize=16)


    flon=np.arange(np.floor(quiv_fact/2).astype('int'),len(FORECAST.lon.values), quiv_fact)
    flat=np.arange(np.floor(quiv_fact/2).astype('int'),len(FORECAST.lat.values), quiv_fact)
    mask=np.meshgrid(flon,flat)

    q1=ax2.quiver(FORECAST.lon.values[flon], 
               FORECAST.lat.values[flat], 
               (np.sin(np.deg2rad(FORECAST.isel(time=time_plot).dir0.values[mask]-180))),
               (np.cos(np.deg2rad(FORECAST.isel(time=time_plot).dir0.values[mask]-180))), 
               FORECAST.isel(time=time_plot).hs0.values[mask], cmap=cmn, clim=(vmin,vmax),
               zorder=3,  transform=ccrs.PlateCarree(), 
#                color='white',
               width=0.003 , headwidth=5, label='Sea'
              )

    q2=ax2.quiver(FORECAST.lon.values[flon], 
               FORECAST.lat.values[flat], 
               (np.sin(np.deg2rad(FORECAST.isel(time=time_plot).dir1.values[mask]-180))),
               (np.cos(np.deg2rad(FORECAST.isel(time=time_plot).dir1.values[mask]-180))), 
               FORECAST.isel(time=time_plot).hs1.values[mask], cmap=cmn, clim=(vmin,vmax),
               zorder=3,  transform=ccrs.PlateCarree(),
#                color='pink',
               width=0.002 , headwidth=5, label='Swell 1'
               )

#     ax2.legend(fontsize=16)
    ax2.set_extent([np.nanmin(FORECAST.lon.values), np.nanmax(FORECAST.lon.values), np.nanmin(FORECAST.lat.values), np.nanmax(FORECAST.lat.values)])
    
    return fig

    
def animate_waves_forecast(FORECAST, coast, sp,anim_interval : int = 300, vmin=[], vmax=[], quiv_fact=1, figsize=[32,13]):
    '''
    Animates waves forecast plot
    '''
    
    if not vmin:
        vmin=np.nanmin(FORECAST.hs.values)
    if not vmax:
        vmax=np.nanmax(FORECAST.hs.values)

    print('generating animation, please wait.')
    
    
    fig = plt.figure(figsize=figsize)
    gs3=gridspec.GridSpec(1,6,hspace=0.1, wspace=0.1)

    time_step=0

    ax1=fig.add_subplot(gs3[0:2],projection='polar')
    z=sp.isel(time=time_step).efth.values
    plot_spectrum(ax1,np.deg2rad(sp.dir.values), sp.freq.values,z, cmap='CMRmap_r', ylim=0.3)
    ax1.set_title('Super Point \n ' + str(sp.isel(time=time_step).time.values.astype('datetime64[h]')), fontsize=16, color='darkmagenta')

    ax2=fig.add_subplot(gs3[2:6],projection = ccrs.PlateCarree(central_longitude=180))
    ax2.add_feature(coast,facecolor='grey', alpha=1, edgecolor='grey', linewidth=1.5, zorder=4)

    cmn = np.vstack((cm.get_cmap('bwr_r', 124)(np.linspace(0.5, 0.96, 20)),
                             cm.get_cmap('jet', 124)(np.linspace(0.2, 0.84, 80)),
                             ))
    cmn = ListedColormap(cmn, name='cmn')

    im=ax2.pcolormesh(FORECAST.lon.values, FORECAST.lat.values, FORECAST.isel(time=time_step).hs.T, vmin=vmin, vmax=vmax, cmap=cmn, transform=ccrs.PlateCarree(),zorder=1, 
                      alpha=0.5, edgecolor='face',
                      antialiased=True, linewidth=0.15)
    
    cbar=ax2.figure.colorbar(im, ax=ax2)
    cbar.set_label('Hs (m)', fontsize=18)
    cbar.ax.tick_params(labelsize=16)


    flon=np.arange(np.floor(quiv_fact/2).astype('int'),len(FORECAST.lon.values), quiv_fact)
    flat=np.arange(np.floor(quiv_fact/2).astype('int'),len(FORECAST.lat.values), quiv_fact)
    mask=np.meshgrid(flon,flat)

    q1=ax2.quiver(FORECAST.lon.values[flon], 
               FORECAST.lat.values[flat], 
               (np.sin(np.deg2rad(FORECAST.isel(time=time_step).dir0.values[mask]-180))),
               (np.cos(np.deg2rad(FORECAST.isel(time=time_step).dir0.values[mask]-180))), 
               FORECAST.isel(time=time_step).hs0.values[mask], cmap=cmn, clim=(vmin,vmax),
               zorder=3,  transform=ccrs.PlateCarree(), 
#                color='white',
               width=0.003 , headwidth=4, label='Sea', edgecolor='navy', linewidth = 0.3
              )

    q2=ax2.quiver(FORECAST.lon.values[flon], 
               FORECAST.lat.values[flat], 
               (np.sin(np.deg2rad(FORECAST.isel(time=time_step).dir1.values[mask]-180))),
               (np.cos(np.deg2rad(FORECAST.isel(time=time_step).dir1.values[mask]-180))), 
               FORECAST.isel(time=time_step).hs1.values[mask], cmap=cmn, clim=(vmin,vmax),
               zorder=3,  transform=ccrs.PlateCarree(),
#                color='pink',
               width=0.002 , headwidth=3, label='Swell 1', edgecolor='navy', linewidth = 0.3
               )
    ax2.set_extent([np.nanmin(FORECAST.lon.values), np.nanmax(FORECAST.lon.values), np.nanmin(FORECAST.lat.values), np.nanmax(FORECAST.lat.values)])
    
    
    def animate(time_step):
        
        ax2.collections=[]
        ax2.patches=[]
        
        #global stp
        #stp.remove()
        
        
        z=sp.isel(time=time_step).efth.values
        plot_spectrum(ax1,np.deg2rad(sp.dir.values), sp.freq.values,z, cmap='CMRmap_r', ylim=0.3)
        
        im=ax2.pcolormesh(FORECAST.lon.values, FORECAST.lat.values, FORECAST.isel(time=time_step).hs.T, vmin=vmin, vmax=vmax, cmap=cmn, transform=ccrs.PlateCarree(),zorder=1, alpha=0.5, edgecolor='face',
                      antialiased=True, linewidth=0.15)
        ax2.add_feature(coast,facecolor='grey', alpha=1, edgecolor='grey', linewidth=1.5, zorder=4)
        
        
#         stp=ax2.streamplot(FORECAST.lon.values, FORECAST.lat.values, (np.sin(np.deg2rad(FORECAST.isel(time=time_step).dir.values-180)).T),
#                            (np.cos(np.deg2rad(FORECAST.isel(time=time_step).dir.values.T-180))), 
#                            zorder=3, linewidth=1.2, color='mediumslateblue',  transform=ccrs.PlateCarree(),
#                            arrowsize=2, density=1.2)
        
        
        q1=ax2.quiver(FORECAST.lon.values[flon], 
               FORECAST.lat.values[flat], 
               (np.sin(np.deg2rad(FORECAST.isel(time=time_step).dir0.values[mask]-180))),
               (np.cos(np.deg2rad(FORECAST.isel(time=time_step).dir0.values[mask]-180))), 
               FORECAST.isel(time=time_step).hs0.values[mask], cmap=cmn, clim=(vmin,vmax),
               zorder=3,  transform=ccrs.PlateCarree(), 
               width=0.003 , headwidth=4, label='Sea', edgecolor='navy', linewidth = 0.3
              )

        q2=ax2.quiver(FORECAST.lon.values[flon], 
               FORECAST.lat.values[flat], 
               (np.sin(np.deg2rad(FORECAST.isel(time=time_step).dir1.values[mask]-180))),
               (np.cos(np.deg2rad(FORECAST.isel(time=time_step).dir1.values[mask]-180))), 
               FORECAST.isel(time=time_step).hs1.values[mask], cmap=cmn, clim=(vmin,vmax),
               zorder=3,  transform=ccrs.PlateCarree(),
               width=0.002 , headwidth=3, label='Swell 1', edgecolor='navy', linewidth = 0.3
               )
        
        tt=ax1.set_title('Super Point \n ' + str(sp.isel(time=time_step).time.values.astype('datetime64[h]')), fontsize=16, color='darkmagenta')


        return [tt]
    
    # Finally, we use the animation module to create the animation
    ani = animation.FuncAnimation(
        fig,
        animate,
        frames = len(FORECAST.time),
        interval = anim_interval,
        repeat = False,
        blit = False,
    )


    plt.close()

    return ani, fig



def plot_resolutions_sl_waves(coast, W_FORECAST, SL_FORECAST, Area_extraction):
    
    figsize=[18,9]

    fig = plt.figure(random.randint(1,50), figsize=figsize)
    gs3=gridspec.GridSpec(1,1)
    ax=fig.add_subplot(gs3[0],projection = ccrs.PlateCarree(central_longitude=180))

    ax.add_feature(coast,facecolor="grey", alpha=0.8, edgecolor='grey', linewidth=1.5, zorder=5)
    xx,yy=np.meshgrid(W_FORECAST.lon.values, W_FORECAST.lat.values)
    ax.plot(xx,yy,'k.', markersize=0.4, transform=ccrs.PlateCarree());
    ax.plot(xx[0],yy[0],'k.', markersize=2, transform=ccrs.PlateCarree(), label='Wave forecast');


    m=ax.plot(SL_FORECAST.lon.values,SL_FORECAST.lat.values,'.', color='coral', markersize=8, transform=ccrs.PlateCarree());
    m=ax.plot(SL_FORECAST.lon.values[0],SL_FORECAST.lat.values[0],'.', color='coral', markersize=8, transform=ccrs.PlateCarree(), label='SL forecast');

#     ax.set_extent([Area_extraction[0], Area_extraction[1], Area_extraction[2], Area_extraction[3]], crs=ccrs.PlateCarree())
    ax.set_extent([W_FORECAST.lon.min(), W_FORECAST.lon.max(), W_FORECAST.lat.min(), W_FORECAST.lat.max(),], crs=ccrs.PlateCarree())
    plt.legend(fontsize=20)
    
    
def animate_components_twl_forecast(W_FORECAST, coast,Area_extraction,  figsize=[29,9], anim_interval : int = 300, vmin=[], vmax=[], quiv_fact=10):
    '''
    Animates waves and sl components of forecast
    '''
    
    if not vmin:
        vmin=np.nanmin(W_FORECAST.hs.values)
    if not vmax:
        vmax=np.nanmax(W_FORECAST.hs.values)

    print('generating animation, please wait.')
    
    time=0
    
    fig = plt.figure(random.randint(1,50), figsize=figsize)
    gs3=gridspec.GridSpec(1,3,hspace=0.1, wspace=0.1)
    
    # SL
    
    ax=fig.add_subplot(gs3[0],projection = ccrs.PlateCarree(central_longitude=180))

    ax.add_feature(coast,facecolor="grey", alpha=0.9, edgecolor='grey', linewidth=1.5, zorder=5)
    
    im=ax.pcolormesh(W_FORECAST.lon.values, W_FORECAST.lat.values, W_FORECAST.isel(time=time).sl.T, vmin=-W_FORECAST.sl.values.max(), vmax=W_FORECAST.sl.values.max(), cmap='RdBu_r', transform=ccrs.PlateCarree(),zorder=1, 
                          alpha=0.9, edgecolor='face',
                          antialiased=True, linewidth=0.15)


    ax.set_extent([Area_extraction[0], Area_extraction[1], Area_extraction[2], Area_extraction[3]], crs=ccrs.PlateCarree())
    plt.colorbar(im, ax=ax,  orientation='horizontal',shrink=0.8, pad=0.05).set_label('SL (m)', fontsize=16)
    
    # Hs
    
    ax1=fig.add_subplot(gs3[1],projection = ccrs.PlateCarree(central_longitude=180))
    ax1.add_feature(coast,facecolor="grey", alpha=0.9, edgecolor='grey', linewidth=1.5, zorder=5)
    cmn = np.vstack((cm.get_cmap('bwr_r', 124)(np.linspace(0.5, 0.96, 20)),
                                 cm.get_cmap('jet', 124)(np.linspace(0.2, 0.84, 80)),
                                 ))
    cmn = ListedColormap(cmn, name='cmn')
    im1=ax1.pcolormesh(W_FORECAST.lon.values, W_FORECAST.lat.values, W_FORECAST.isel(time=time).hs.T, vmin=0, vmax=W_FORECAST.hs.values.max(), cmap=cmn, transform=ccrs.PlateCarree(),zorder=1)

    flon=np.arange(np.floor(quiv_fact/2).astype('int'),len(W_FORECAST.lon.values), quiv_fact)
    flat=np.arange(np.floor(quiv_fact/2).astype('int'),len(W_FORECAST.lat.values), quiv_fact)
    mask=np.meshgrid(flon,flat)
    
    q1=ax1.quiver(W_FORECAST.lon.values[flon], 
               W_FORECAST.lat.values[flat], 
               (np.sin(np.deg2rad(W_FORECAST.isel(time=time).dir.values[mask]-180))),
               (np.cos(np.deg2rad(W_FORECAST.isel(time=time).dir.values[mask]-180))), 
               color='grey',
               zorder=3,  transform=ccrs.PlateCarree(), 
               width=0.003 , headwidth=5, label='Sea'
              )
    
    
    ax1.set_extent([Area_extraction[0], Area_extraction[1], Area_extraction[2], Area_extraction[3]], crs=ccrs.PlateCarree())
    plt.colorbar(im1, ax=ax1,  orientation='horizontal',shrink=0.8, pad=0.05).set_label('Hs (m)', fontsize=16)

    ## TWL

    ax2=fig.add_subplot(gs3[2],projection = ccrs.PlateCarree(central_longitude=180))

    ax2.add_feature(coast,facecolor="grey", alpha=0.9, edgecolor='grey', linewidth=1.5, zorder=5)

    im2=ax2.pcolormesh(W_FORECAST.lon.values, W_FORECAST.lat.values, W_FORECAST.isel(time=time).twl.T, vmin=0, vmax=W_FORECAST.twl.values.max(), cmap='rainbow', transform=ccrs.PlateCarree(),zorder=1, 
                          alpha=1, edgecolor='face',
                          antialiased=True, linewidth=0.15)


    ax2.set_extent([Area_extraction[0], Area_extraction[1], Area_extraction[2], Area_extraction[3]], crs=ccrs.PlateCarree())
    plt.colorbar(im2, ax=ax2,  orientation='horizontal',shrink=0.8, pad=0.05).set_label('TWL (m)', fontsize=16)
    
    tt=ax1.set_title('Time: ' + str(W_FORECAST.isel(time=time).time.values.astype('datetime64[s]')), fontsize=22, color='darkmagenta', pad=50)

    
    def animate(time):
        
        ax1.collections=[]
        ax1.patches=[]
          
        im=ax.pcolormesh(W_FORECAST.lon.values, W_FORECAST.lat.values, W_FORECAST.isel(time=time).sl.T, vmin=-W_FORECAST.sl.values.max(), vmax=W_FORECAST.sl.values.max(), cmap='RdBu_r', transform=ccrs.PlateCarree(),zorder=1, alpha=0.9)
        
        
        im1=ax1.pcolormesh(W_FORECAST.lon.values, W_FORECAST.lat.values, W_FORECAST.isel(time=time).hs.T, vmin=0, vmax=W_FORECAST.hs.values.max(), cmap=cmn, transform=ccrs.PlateCarree(),zorder=1)
        q1=ax1.quiver(W_FORECAST.lon.values[flon], 
               W_FORECAST.lat.values[flat], 
               (np.sin(np.deg2rad(W_FORECAST.isel(time=time).dir.values[mask]-180))),
               (np.cos(np.deg2rad(W_FORECAST.isel(time=time).dir.values[mask]-180))), 
               color='grey',
               zorder=3,  transform=ccrs.PlateCarree(), 
               width=0.003 , headwidth=5, label='Sea'
              )
        
        im2=ax2.pcolormesh(W_FORECAST.lon.values, W_FORECAST.lat.values, W_FORECAST.isel(time=time).twl.T, vmin=0, vmax=W_FORECAST.twl.values.max(), cmap='rainbow', transform=ccrs.PlateCarree(),zorder=1)
        
        tt=ax1.set_title('Time: ' + str(W_FORECAST.isel(time=time).time.values.astype('datetime64[s]')), fontsize=22, color='darkmagenta', pad=50)

        return [tt]

      

    # Finally, we use the animation module to create the animation
    ani = animation.FuncAnimation(
        fig,
        animate,
        frames = len(W_FORECAST.time),
        interval = anim_interval,
        repeat = True,
        blit = False,
    )


    plt.close()

    return ani, fig



def animate_twl_forecast(W_FORECAST, coast,Area_extraction,  anim_interval : int = 300, vmin=[], vmax=[], figsize=[15,9]):
    '''
    Animates TWL forecast
    '''
    
    if not vmin:
        vmin=np.nanmin(W_FORECAST.hs.values)
    if not vmax:
        vmax=np.nanmax(W_FORECAST.hs.values)

    print('generating animation, please wait.')
    
    time=0
    
    figsize=figsize

    x_grid, y_grid=np.meshgrid(W_FORECAST.lon.values, W_FORECAST.lat.values)

    fig = plt.figure(random.randint(1,50), figsize=figsize)
    
    gs3=gridspec.GridSpec(1,1,hspace=0.1, wspace=0.1)
    ax=fig.add_subplot(gs3[0],projection = ccrs.PlateCarree(central_longitude=180))


    ax.add_feature(coast,facecolor="grey", alpha=0.9, edgecolor='grey', linewidth=1.5, zorder=5)
    # im=ax.scatter(x_grid.T, y_grid.T, 10, W_FORECAST.isel(time=time).sl.values,vmin=-1, vmax=1, cmap='RdBu_r', transform=ccrs.PlateCarree())

    ax.add_feature(coast,facecolor="grey", alpha=1, edgecolor='grey', linewidth=1.5, zorder=5)
    levels=np.arange(W_FORECAST.twl.values.min()-0.1,W_FORECAST.twl.values.max()+0.1,0.05)    
    im=ax.contourf(W_FORECAST.lon.values, W_FORECAST.lat.values, W_FORECAST.isel(time=time).twl.T, levels, vmin=W_FORECAST.twl.values.min(), vmax=W_FORECAST.twl.values.max(), cmap='rainbow', transform=ccrs.PlateCarree(),zorder=1, alpha=1)

    ax.set_extent([Area_extraction[0], Area_extraction[1], Area_extraction[2], Area_extraction[3]], crs=ccrs.PlateCarree())
    plt.colorbar(im, ax=ax,  orientation='vertical',shrink=0.8, pad=0.05).set_label('TWL (m)', fontsize=16)

    tt=ax.set_title('Time: ' + str(W_FORECAST.isel(time=time).time.values.astype('datetime64[s]')), fontsize=18, color='darkmagenta')

    
    def animate(time):
        
        ax.collections=[]
        ax.patches=[]
                
        im=ax.contourf(W_FORECAST.lon.values, W_FORECAST.lat.values, W_FORECAST.isel(time=time).twl.T, levels, vmin=W_FORECAST.twl.values.min(), vmax=W_FORECAST.twl.values.max(), cmap='rainbow', transform=ccrs.PlateCarree(),zorder=1, alpha=1)
        
        tt=ax.set_title('Time: ' + str(W_FORECAST.isel(time=time).time.values.astype('datetime64[s]')), fontsize=18, color='darkmagenta', pad=30)

        return [tt]

    # Finally, we use the animation module to create the animation
    ani = animation.FuncAnimation(
        fig,
        animate,
        frames = len(W_FORECAST.time),
        interval = anim_interval,
        repeat = True,
        blit = False,
    )

    plt.close()

    return ani, fig


def plot_hindcast_time(time_stats, coast, h_min=0, h_max=4, t_min=6, t_max=18):
    
    fig = plt.figure(figsize=[22,23])
    gs3=gridspec.GridSpec(2,2,hspace=0.01, wspace=0.01) 

    var_s='Hs'
    ax=fig.add_subplot(gs3[0],projection = ccrs.PlateCarree(central_longitude=180))
    Plot_hindcast_maps(time_stats, coast, var_s, p=[], vmin=h_min, vmax=h_max, ax=ax)
    ax.set_title(var_s, fontsize=16, color='firebrick', fontweight='bold')

    var_s='Tp'
    ax=fig.add_subplot(gs3[1],projection = ccrs.PlateCarree(central_longitude=180))
    Plot_hindcast_maps(time_stats, coast, var_s, p=[], vmin=t_min, vmax=t_max, ax=ax)
    ax.set_title(var_s, fontsize=16, color='firebrick', fontweight='bold')

    
    
def Plot_Diferences_Real_Cluster(sp_regular,predict_cluster_sorted, efth_cluster, pos, ylim=0.5, vmax=0.7, ymax=0.3, cmap='magma_r'):
    
    fig = plt.figure(figsize=[18,7])


    gs3=gridspec.GridSpec(1,3,hspace=0.01, wspace=0.01)
    ax0=fig.add_subplot(gs3[0],projection='polar')

    x=np.deg2rad(sp_regular.dir.values)
    y=sp_regular.freq.values
    z= np.sqrt(sp_regular.efth.values[pos,:,:])
    p2=plot_spectrum(ax0,x,y, z, vmin=0, vmax=vmax, ylim=np.nanmax(sp_regular.freq), remove_axis=1, cmap=cmap)
    ax0.set_ylim([0, ylim])
    ax0.set_title('Real ( Hs = '+ str(np.round(sp_regular.efth.spec.hs().values[pos], 3)) + ')', fontsize=18)

    ax1=fig.add_subplot(gs3[1],projection='polar')
    z1= np.sqrt(efth_cluster.efth.values[np.int64(predict_cluster_sorted[pos]),:,:])
    p3=plot_spectrum(ax1,x,y, z1, vmin=0, vmax=vmax, ylim=np.nanmax(sp_regular.freq), remove_axis=1, cmap=cmap)
    ax1.set_ylim([0, ylim])
    ax1.set_title('Cluster ( Hs = '+ str(np.round(efth_cluster.efth.spec.hs().values[np.int64(predict_cluster_sorted[pos])], 3)) + ')', fontsize=18)

    ax2=fig.add_subplot(gs3[2],projection='polar')
    z1= np.sqrt(efth_cluster.efth.values[np.int64(predict_cluster_sorted[pos]),:,:])
    p4=plot_spectrum(ax2,x,y, z1-z, vmin=-0.2, vmax=0.2, ylim=np.nanmax(sp_regular.freq), remove_axis=1, cmap='RdBu_r')
    ax1.set_ylim([0, ylim])
    ax2.set_title('Cluster - Real', fontsize=18)
    gs3.tight_layout(fig, rect=[[], [], 0.85, []])

    gs4=gridspec.GridSpec(1,2)
    ax0=fig.add_subplot(gs4[0])
    plt.colorbar(p2,cax=ax0, shrink=0.7)
    ax0.set_ylabel('Sqrt (Efth)', fontsize=14)

    ax0=fig.add_subplot(gs4[1])
    plt.colorbar(p4,cax=ax0, shrink=0.7)
    ax0.set_ylabel('Clust-Real', fontsize=14)

    gs4.tight_layout(fig, rect=[0.85, 0.2, 1, 0.8])