# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 00:30:50 2022

@author: moritzw
"""
import pandas as pd
df = pd.read_csv('Vaitupu_Oceanside_Profilesv2.csv')

bad_transect = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,50,51,208,247,248,249,250,251,252,253,254,255,256,257,258,258,259,260,261]

for i in range(len(bad_transect)):
    bad_ix = df.LINE_ID[df.Transect == bad_transect[i]]
    df = df.drop(bad_ix.index)

df.to_csv('Vaitupu_Oceanside_Profilesv3.csv',index = False)
