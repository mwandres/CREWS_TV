from rbf import RBF_Reconstruction, RBF_Validation, fminbound, CostEps, CalcRBF_Coeff
import numpy as np
from mda import Normalize
import pandas as pd
import time

def get_rbf_coefficients_and_save_as_pkl(training_dataset,
        ix_scalar_subset, ix_directional_subset,
        response_dataset,
        ix_scalar_target, out_name):
    # normalize subset and dataset
    # dataset_norm, mins, maxs = Normalize(
    #     target_data, ix_scalar_subset, ix_directional_subset)
    
    subset_norm,mins, maxs = Normalize(
        training_dataset, ix_scalar_subset, ix_directional_subset)
    sigma_min=0.1
    sigma_max=0.7
    df = pd.DataFrame(data = {'rbf_coeff' : [], 'opt_sigma' : [], 'mins' : [], 'maxs' : []})
    df.append({'mins' : mins, 'maxs' : maxs},ignore_index = True)
    for ix in ix_scalar_target:
        v = response_dataset[:,ix]
        t0 = time.time()
        # minimize RBF cost function
        opt_sigma = fminbound(
            CostEps, sigma_min, sigma_max, args=(subset_norm.T, v)
        )
        t1 = time.time()  # optimization time
        rbf_coeff, _ = CalcRBF_Coeff(opt_sigma, subset_norm.T, v)
        df = df.append(
            {'rbf_coeff' : rbf_coeff, 
            'opt_sigma' : opt_sigma },ignore_index = True
            )
        status = 'Getting RBF coefficients for ' + str(ix) + '/' + str(len(ix_scalar_target))
        print(status)
        print(
            'ix_scalar: {0},  optimization: {1:.2f} |, sigma_opt: {2:.2f}'.format(
                ix, t1-t0, opt_sigma)
            )        
    df.to_pickle(out_name)
    return(df)


import os
os.chdir('F:/XBeach/Tuvalu/1D_sims/Nui_forecast/')
transects_csv = 'F:/XBeach/Tuvalu/1D_sims/Nui_forecast/Nui/Nui_Oceanside_Profilesv4.csv'
transects = pd.read_csv(transects_csv,delimiter=',',skiprows=0)
transect_labels =pd.unique(transects.LINE_ID)
transect_labels = transect_labels.astype(int)

for i in range(len(transect_labels)):
    #print(transect_labels[i])
 
#i = 0

    in_file_name = 't_' + str(transect_labels[i]) + '_results.csv'
    out_name = 't_' + str(transect_labels[i]) + '_rbf_coefficients_max_TWL_and_inun_ext.pkl'
    print(in_file_name)
    # out_name_2 = 't_' + str(transect_labels[i]) + '_rbf_coefficients_transect.pkl'
    
    # load transect data
    data = np.loadtxt(in_file_name,delimiter=',',skiprows=0)
    training_dataset = data[:,0:4]
    response_dataset_1 = data[:,4:6]
    response_dataset_2 = data[:,6:]
    
    # subset - scalar / directional indexes
    ix_scalar_subset = [0,1,3]      # scalar (hs, tp)
    ix_directional_subset = [2]   # directional (dir)
    
    # target - scalar / directional indexes
    ix_scalar_target = [0,1]      # scalar (Hsig, Tpsmoo)
    ix_directional_target = []   # directional (Dir)
    ix_scalar_target_2 = list(np.arange(0,len(response_dataset_2[0,:])))
    ix_directional_target_2 = []
     
    #calculate RBF coefficients and save as pickle
    get_rbf_coefficients_and_save_as_pkl(
        training_dataset,
        ix_scalar_subset, ix_directional_subset,
        response_dataset_1,
        ix_scalar_target, out_name
        )
    print(out_name)

 # get_rbf_coefficients_and_save_as_pkl(
 #     training_dataset,
 #     ix_scalar_subset, ix_directional_subset,
 #     response_dataset_2,
 #     ix_scalar_target_2, out_name_2
 #     )


# output[:, ix] = RBF_Interpolation(
#     opt_sigma, rbf_coeff, subset_norm.T, dataset_norm.T)



# # reconstruct transect
# output_2_inundation_transect = RBF_Reconstruction(
#     training_dataset, ix_scalar_subset, ix_directional_subset,
#     response_dataset_2, ix_scalar_target_2, ix_directional_target_2,
#     target_data,
# )

#reconstruct TWL nearshore and inundation extent
# output_1 = RBF_Reconstruction(
#     training_dataset, ix_scalar_subset, ix_directional_subset,
#     response_dataset_1, ix_scalar_target, ix_directional_target,
#     target_data,
# )

# output_1 = RBF_Validation(
#     training_dataset, ix_scalar_subset, ix_directional_subset,
#     response_dataset_1, ix_scalar_target, ix_directional_target,
# )
# print(output_1)
