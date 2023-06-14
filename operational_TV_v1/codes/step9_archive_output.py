# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 14:10:27 2021

@author: judithg
"""

# Archiving the output.mat

import os, shutil
import glob

def archive_output(now):
# function to save the swan output.mat file to the archive folder
    #floutname='../archives/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") + now.strftime("%H")  +'_output.mat'
    #fout_name = '../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") + now.strftime("%H")  + '/output.mat'
    #shutil.copy(fout_name, floutname)
    #Archieve inundation forecast
    dirname = '../archives/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") + now.strftime("%H")
    if os.path.exists(dirname):
        print('exists')
    else:
        os.mkdir(dirname,775)
    inun_dir='../inundation/Figures'
    #shutil.copy(inun_dir, dirname)
    for file in glob.glob(os.path.join(inun_dir,"*.png")):
        shutil.copy2(file,dirname)
    for file in glob.glob(os.path.join('../inundation/Flood_risk',"*.csv")):
        shutil.copy2(file,dirname)
    
    floutname='../archives/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") + now.strftime("%H")  +'/output.mat'
    fout_name = '../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") + now.strftime("%H")  + '/output.mat'
    shutil.copy(fout_name, floutname)
    print('done')
    #ends
    return()
    
