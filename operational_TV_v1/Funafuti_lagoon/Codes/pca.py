# common
from datetime import datetime

# pip
import numpy as np
import xarray as xr
from sklearn.decomposition import PCA


def spatial_gradient(xdset, var_name):
    '''
    Calculate spatial gradient

    xdset:
        (longitude, latitude, time), var_name

    returns xdset with new variable "var_name_gradient"
    '''

    # TODO:check/ ADD ONE ROW/COL EACH SIDE
    var_grad = np.zeros(xdset[var_name].shape)

    lat = xdset.latitude.values

    for it in range(len(xdset.time)):
        var_val = xdset[var_name].isel(time=it).values

        # calculate gradient (matrix)
        m_c = var_val[1:-1,1:-1]
        m_l = np.roll(var_val, -1, axis=1)[1:-1,1:-1]
        m_r = np.roll(var_val, +1, axis=1)[1:-1,1:-1]
        m_u = np.roll(var_val, -1, axis=0)[1:-1,1:-1]
        m_d = np.roll(var_val, +1, axis=0)[1:-1,1:-1]
        m_phi = np.pi*np.abs(lat)/180.0
        m_phi = m_phi[1:-1]

        dpx1 = (m_c - m_l)/np.cos(m_phi[:,None])
        dpx2 = (m_r - m_c)/np.cos(m_phi[:,None])
        dpy1 = m_c - m_d
        dpy2 = m_u - m_c

        vg = (dpx1**2+dpx2**2)/2 + (dpy1**2+dpy2**2)/2
        var_grad[it, 1:-1, 1:-1] = vg

    # store gradient
    xdset['{0}_gradient'.format(var_name)]= (
        ('time', 'latitude', 'longitude'), var_grad)

    return xdset
    
    
def dynamic_estela_predictor(xdset, var_name, estela_D, n_days=None):
    '''
    Generate dynamic predictor using estela

    xdset:
        (time, latitude, longitude), var_name, mask

    returns similar xarray.Dataset with variables:
        (time, latitude, longitude), var_name_comp
        (time, latitude, longitude), var_name_gradient_comp
    '''

    # first day is estela max
    if n_days:
    	first_day = n_days
    else:
    	first_day = int(np.floor(np.nanmax(estela_D)))+1

    # output will start at time=first_day
    shp = xdset[var_name].shape
    comp_shape = (shp[0]-first_day, shp[1], shp[2])
    var_comp = np.ones(comp_shape) * np.nan
    var_grd_comp = np.ones(comp_shape) * np.nan

    # get data using estela for each cell
    for i_lat in range(len(xdset.latitude)):
        for i_lon in range(len(xdset.longitude)):
            ed = estela_D[i_lat, i_lon]
            if not np.isnan(ed):

                # mount estela displaced time array 
                i_times = np.arange(
                    first_day, len(xdset.time)
                ) - np.int(ed)

                # select data from displaced time array positions
                xdselec = xdset.isel(
                    time = i_times,
                    latitude = i_lat,
                    longitude = i_lon)

                # get estela predictor values
                var_comp[:, i_lat, i_lon] = xdselec[var_name].values
                var_grd_comp[:, i_lat, i_lon] = xdselec['{0}_gradient'.format(var_name)].values

    # return generated estela predictor
    return xr.Dataset(
        {
            '{0}_comp'.format(var_name):(
                ('time','latitude','longitude'), var_comp),
            '{0}_gradient_comp'.format(var_name):(
                ('time','latitude','longitude'), var_grd_comp),

        },
        coords = {
            'time':xdset.time.values[first_day:],
            'latitude':xdset.latitude.values,
            'longitude':xdset.longitude.values,
        }
    )


# def PCA_EstelaPred(xds, pred_name):
#     '''
#     Principal component analysis
#     method: custom for estela predictor

#     xds:
#         (time, latitude, longitude), pred_name_comp | pred_name_gradient_comp

#     returns a xarray.Dataset containing PCA data: PCs, EOFs, variance
#     '''

#     # estela predictor and estela gradient predictor
#     try:
#         pred_est_var = xds['{0}_comp'.format(pred_name)]
#         pred_est_grad = xds['{0}_gradient_comp'.format(pred_name)]
#     except:
#         pred_est_var = xds['{0}'.format(pred_name)]
#         pred_est_grad = xds['{0}_gradient'.format(pred_name)]

#     # use data inside timeframe
#     dp_var = pred_est_var.values
#     dp_grd = pred_est_grad.values
#     shape_grid = dp_var[0].shape  # needed to handle data after PCs

#     # unravel and join var and grad data 
#     dp_ur = np.nan * np.ones(
#         (dp_var.shape[0], 2*dp_var.shape[1]*dp_var.shape[2])
#     )

#     # we use .T to equal matlab
#     for ti in range(dp_ur.shape[0]):
#         dp_ur[ti,:] = np.concatenate(
#             [np.ravel(dp_var[ti].T) , np.ravel(dp_grd[ti].T)]
#         )

#     # remove nans from predictor    
#     data_pos = ~np.isnan(dp_ur[0,:])
#     clean_row = dp_ur[0, data_pos]
#     dp_ur_nonan = np.nan * np.ones(
#         (dp_ur.shape[0], len(clean_row))
#     )
#     for ti in range(dp_ur.shape[0]):
#         dp_ur_nonan[ti,:] = dp_ur[ti, data_pos]

#     # standarize predictor
#     pred_mean = np.mean(dp_ur_nonan, axis=0)
#     pred_std = np.std(dp_ur_nonan, axis=0)
#     pred_norm = (dp_ur_nonan[:,:] - pred_mean) / pred_std
#     pred_norm[np.isnan(pred_norm)] = 0

#     # principal components analysis
#     ipca = PCA(n_components=min(pred_norm.shape[0], pred_norm.shape[1]))
#     PCs = ipca.fit_transform(pred_norm)

#     # return dataset
#     a = xr.Dataset(
#         {
#             'PCs': (('time', 'n_components'), PCs),
#             'EOFs': (('n_components','n_features'), ipca.components_),
#             'variance': (('n_components',), ipca.explained_variance_),

#             'pred_mean': (('n_features',), pred_mean),
#             'pred_std': (('n_features',), pred_std),

#             'pred_lon': (('n_lon',), xds.longitude.values[:]),
#             'pred_lat': (('n_lat',), xds.latitude.values[:]),
#             'pred_time': (('time',), xds.time.values[:]),
#             'pred_data_pos':(('n_points',), data_pos)
#         },

#         attrs = {
#             'method': 'gradient + estela',
#             'pred_name': pred_name,
#         }
#     )

#     return a, ipca

def standardise_predictor_prep(xds, pred_name):
    '''
    Predictor dataset is standardized for PCA
    '''

    # estela predictor and estela gradient predictor
    try:
        pred_est_var = xds['{0}_comp'.format(pred_name)]
        pred_est_grad = xds['{0}_gradient_comp'.format(pred_name)]
    except:
        pred_est_var = xds['{0}'.format(pred_name)]
        pred_est_grad = xds['{0}_gradient'.format(pred_name)]

    # use data inside timeframe
    dp_var = pred_est_var.values
    dp_grd = pred_est_grad.values
    shape_grid = dp_var[0].shape  # needed to handle data after PCs
    
    # unravel and join var and grad data 
    dp_ur = np.nan * np.ones(
        (dp_var.shape[0], 2*dp_var.shape[1]*dp_var.shape[2])
    )

    # we use .T to equal matlab
    for ti in range(dp_ur.shape[0]):
        dp_ur[ti,:] = np.concatenate(
            [np.ravel(dp_var[ti].T) , np.ravel(dp_grd[ti].T)]
        )

    # remove nans from predictor    
    data_pos = ~np.isnan(dp_ur[0,:])
    clean_row = dp_ur[0, data_pos]
    dp_ur_nonan = np.nan * np.ones(
        (dp_ur.shape[0], len(clean_row))
    )
    for ti in range(dp_ur.shape[0]):
        dp_ur_nonan[ti,:] = dp_ur[ti, data_pos]

    # standarize predictor
    pred_mean = np.mean(dp_ur_nonan, axis=0)
    pred_std = np.std(dp_ur_nonan, axis=0)
    pred_norm = (dp_ur_nonan[:,:] - pred_mean) / pred_std
    pred_norm[np.isnan(pred_norm)] = 0
    

    return dp_ur_nonan, data_pos

def standardise_predictor(xds, pred_name_ls):
    '''
    Predictor dataset is standardized for PCA
    '''
   
    dp_ur_nonan, data_pos = standardise_predictor_prep(xds, pred_name_ls)
    
    # standarize predictor
    pred_mean = np.mean(dp_ur_nonan, axis=0)
    pred_std = np.std(dp_ur_nonan, axis=0)
    pred_norm = (dp_ur_nonan[:,:] - pred_mean) / pred_std
    pred_norm[np.isnan(pred_norm)] = 0      # RuntimeWarning: invalid value encountered in true_divide !!!!!????
   
    # return dataset
    xds_norm =  xr.Dataset(
        {  
            'pred_norm': (('time', 'n_components'), pred_norm),
            'pred_mean': (('n_features',), pred_mean),
            'pred_std': (('n_features',), pred_std),
        },

        coords = {
            'time': xds.time.values[:],
            'lon': xds.longitude.values[:],
            'lat': xds.latitude.values[:]
        },
         attrs = {
            'method': 'gradient + estela + sea mask',
            'pred_name': pred_name_ls,
        }
    )
        
    return data_pos, xds_norm


def spatial_gradient(xdset, var_name):
    '''
    Calculate spatial gradient

    xdset:
        (longitude, latitude, time), var_name

    returns xdset with new variable "var_name_gradient"
    '''


    # TODO:check/ ADD ONE ROW/COL EACH SIDE
    var_grad = np.zeros(xdset[var_name].shape)

    lat = xdset.latitude.values

    for it in range(len(xdset.time)):
        var_val = xdset[var_name].isel(time=it).values

        # calculate gradient (matrix)
        m_c = var_val[1:-1,1:-1]
        m_l = np.roll(var_val, -1, axis=1)[1:-1,1:-1]
        m_r = np.roll(var_val, +1, axis=1)[1:-1,1:-1]
        m_u = np.roll(var_val, -1, axis=0)[1:-1,1:-1]
        m_d = np.roll(var_val, +1, axis=0)[1:-1,1:-1]
        m_phi = np.pi*np.abs(lat)/180.0
        m_phi = m_phi[1:-1] 

        dpx1 = (m_c - m_l)/np.cos(m_phi[:,None])
        dpx2 = (m_r - m_c)/np.cos(m_phi[:,None])
        dpy1 = m_c - m_d
        dpy2 = m_u - m_c

        vg = (dpx1**2+dpx2**2)/2 + (dpy1**2+dpy2**2)/2
        var_grad[it, 1:-1, 1:-1] = vg

    # store gradient
    xdset['{0}_gradient'.format(var_name)]= (
        ('time', 'latitude', 'longitude'), var_grad)

    return xdset

def PCA_EstelaPred(xds, pred_name_ls):
    '''
    Principal component analysis
    method: custom for estela predictor

    xds:
        (time, latitude, longitude), pred_name | pred_name_gradient

    returns a xarray.Dataset containing PCA data: PCs, EOFs, variance
    '''

    # standardise predictor data
    data_pos, xds_norm = standardise_predictor(xds, pred_name_ls)
   
    # normalized predictor
    pred_norm = xds_norm.pred_norm.values
   
    # principal components analysis
    ipca = PCA(n_components=min(pred_norm.shape[0], pred_norm.shape[1]))
    PCs = ipca.fit_transform(pred_norm)

    # return data
    
    a = xr.Dataset(
        {
            'PCs': (('time', 'n_components'), PCs),
            'EOFs': (('n_components','n_features'), ipca.components_),
            'variance': (('n_components',), ipca.explained_variance_),

            'pred_mean': (('n_features',), xds_norm.pred_mean),
            'pred_std': (('n_features',), xds_norm.pred_std),

            'pred_lon': (('n_lon',), xds.longitude.values[:]),
            'pred_lat': (('n_lat',), xds.latitude.values[:]),
            'pred_time': (('time',), xds.time.values[:]),
            'pred_data_pos':(('n_points',), data_pos)
        },

        attrs = {
            'method': 'gradient + estela',
            'pred_name': xds_norm.pred_name,
        }
    )

    return a, ipca