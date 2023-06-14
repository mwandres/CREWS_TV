# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 00:30:50 2022

@author: moritzw
"""
import pandas as pd
df = pd.read_csv('Nui_Oceanside_Profilesv3.csv')

bad_transect = [0,80,82,83,84,85,86,88,89,90,91,92,119,125,142,143,144,145,
                146,147,148,149,150,151,152,153,154,155,156,157,158,165,
                166,167,178,179,180,181,182,183,184,185,186,213,214,215,
                216,217,218,219,220,221,222,223,224,225,226,227,228,229,
                230,231,232,233,235,237]

for i in range(len(bad_transect)):
    bad_ix = df.LINE_ID[df.LINE_ID == bad_transect[i]]
    df = df.drop(bad_ix.index)

df.to_csv('Nui_Oceanside_Profilesv4.csv',index = False)
