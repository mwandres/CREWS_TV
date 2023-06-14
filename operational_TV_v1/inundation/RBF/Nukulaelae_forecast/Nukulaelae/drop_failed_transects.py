# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 00:30:50 2022

@author: moritzw
"""
import pandas as pd
df = pd.read_csv('Nukulaelae_Oceanside_Profilesv3.csv')

bad_transect = [5,12,13,51,101,102,118,119,120,121,122,240,241,242,243,262,263,264,266,273,274,275,276,277,278,279,294,295,307,309,310,311,312,313,314,315]

for i in range(len(bad_transect)):
    bad_ix = df.LINE_ID[df.LINE_ID == bad_transect[i]]
    df = df.drop(bad_ix.index)

df.to_csv('Nukulaelae_Oceanside_Profilesv4.csv',index = False)
