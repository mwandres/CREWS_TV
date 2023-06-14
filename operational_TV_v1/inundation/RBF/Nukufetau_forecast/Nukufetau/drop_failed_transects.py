# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 00:30:50 2022

@author: moritzw
"""
import pandas as pd
df = pd.read_csv('Nukufetau_Oceanside_Profilesv3.csv')

bad_transect = [436]

for i in range(len(bad_transect)):
    bad_ix = df.LINE_ID[df.LINE_ID == bad_transect[i]]
    df = df.drop(bad_ix.index)

df.to_csv('Nukufetau_Oceanside_Profilesv4.csv',index = False)
