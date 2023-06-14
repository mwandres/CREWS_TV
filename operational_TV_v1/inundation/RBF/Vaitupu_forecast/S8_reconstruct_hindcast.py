# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 04:05:12 2021

@author: moritzw
"""
import scipy.io
import numpy as np
import pandas as pd
from rbf import RBF_Interpolation, Normalize
from matplotlib import pyplot as plt
import os

def calculate_thresholds_from_historic_data(hist_data_path,training_data_path,rbf_coeff_path,output_thresholds_path):
    #change period factor when looping through transects
    mat = scipy.io.loadmat(hist_data_path)
    
    # build a list of keys and values for each entry in the structure
    vals = mat['data']['hist'][0,0] #<-- set the array you want to access. 
    keys = mat['data']['hist'][0,0].dtype.descr
    
    # Assemble the keys and values into variables with the same name as that used in MATLAB
    # for i in range(len(keys)):
    #     key = keys[i][0]
    #     val = np.squeeze(vals[key][0][0])  # squeeze is used to covert matlat (1,n) arrays into numpy (1,) arrays. 
    #     exec(key + '=val')
    Hs = np.squeeze(vals['Hs'][0][0])    
    Per = np.squeeze(vals['Per'][0][0])
    Dir = np.squeeze(vals['Dir'][0][0])    
    MSLA = np.squeeze(vals['MSLA'][0][0])
    
    period_fac = mat['data']['period_fac'][0,0]
    Tp = Per * period_fac[0,0] # change factor when calculating different transect
    Hs = Hs.reshape((len(Hs),1))
    Tp = Tp.reshape((len(Tp),1))
    Dir = Dir.reshape((len(Dir),1))
    MSLA = MSLA.reshape((len(MSLA),1))    
    
    target_data = np.stack((Hs,Tp,Dir,MSLA),axis=1).reshape((len(Hs),4))
    
    # data = np.loadtxt('t_1_results.csv',delimiter=',',skiprows=0)
    # target_data = data[:,0:4]
    
    ##
    # load training data
    training_data_csv = np.loadtxt(training_data_path,delimiter=',',skiprows=0)
    training_dataset = training_data_csv[:,0:4]
    # subset - scalar / directional indexes
    ix_scalar = [0,1,3]      # scalar (hs, tp)
    ix_directional = [2]   # directional (dir)
    ##
    df = pd.read_pickle(rbf_coeff_path)

    #target_data =np.ndarray(shape=(1,4), buffer=np.array([3.,10.,150.,1.]), dtype=float, order='F')
    
    # normalize subset and dataset
    subset_norm, mins, maxs = Normalize(
        training_dataset, ix_scalar, ix_directional)
    
    dataset_norm, _, _ = Normalize(
        target_data, ix_scalar, ix_directional,mins,maxs)
    
    opt_sigma = df.opt_sigma[0]
    rbf_coeff = df.rbf_coeff[0]
    output_twl_nearshore = RBF_Interpolation(opt_sigma, rbf_coeff, subset_norm.T, dataset_norm.T)
    output_twl_nearshore = np.nan_to_num(output_twl_nearshore)
    
    
    output_twl_nearshore_sorted = np.sort(output_twl_nearshore)
    Promedioz = np.arange(1,len(output_twl_nearshore)+1)/(len(output_twl_nearshore)+1)
    
    Tamedioz = 1/(1-Promedioz)
    lbda = len(output_twl_nearshore)/40
    
    ari_x = Tamedioz/lbda
    
    a = np.array((ari_x[-80:],output_twl_nearshore_sorted[-80:])).T
    np.savetxt(output_thresholds_path, a, delimiter=",")
    
    # plt.figure()
    # plt.semilogx(ari_x,output_twl_nearshore_sorted)
    return(a)


os.chdir('G:/XBeach/Tuvalu/1D_sims/Vaitupu_forecast/Vaitupu_forecast/')

transects_csv = 'G:/XBeach/Tuvalu/1D_sims/Vaitupu_forecast/Vaitupu_forecast/Vaitupu/Vaitupu_Oceanside_Profilesv3.csv'
transects = pd.read_csv(transects_csv,delimiter=',',skiprows=0)
transect_labels =pd.unique(transects.Transect)
transect_labels = transect_labels.astype(int)

for i in range(len(transect_labels)):
    hist_data_path = 'F:/XBeach/Tuvalu/1D_sims/Vaitupu/t_'+str(transect_labels[i])+'_base/hist_and_ext_data.mat'
    training_data_path = 'G:/XBeach/Tuvalu/1D_sims/Vaitupu_forecast/Vaitupu_forecast/t_'+str(transect_labels[i])+'_results.csv'
    rbf_coeff_path = 'G:/XBeach/Tuvalu/1D_sims/Vaitupu_forecast/Vaitupu_forecast/t_'+str(transect_labels[i])+'_rbf_coefficients_max_TWL_and_inun_ext.pkl'
    output_thresholds_path = 'G:/XBeach/Tuvalu/1D_sims/Vaitupu_forecast/Vaitupu_forecast/t_'+str(transect_labels[i])+'_thresholds.csv'

    a = calculate_thresholds_from_historic_data(hist_data_path,training_data_path,rbf_coeff_path,output_thresholds_path)
    print(output_thresholds_path)
    print(str(a[0,0]) + ' year ARI is ' + str(a[0,1]))
    print(str(a[72,0]) + ' year ARI is ' + str(a[72,1]))


