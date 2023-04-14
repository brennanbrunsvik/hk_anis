#%%
import numpy as np 
import matplotlib.pyplot as plt 
import obspy 
import pandas as pd 
import csv


#%%
fname_sum  = 'ears_sta_summary.csv'
df = pd.read_csv(fname_sum)
df = df[(df.Long > -125) * (df.Long < -65) * (df.Lat > 25) * (df.Lat < 50)]; 

thickstr = df['Est.Thick'].to_numpy()
stdstr = df['StdDev'].to_numpy()

thick = []
std = []
for ista in range(len(thickstr)): 
    mohoi = thickstr[ista]
    mohoi = float(mohoi[0:2]) 
    thick.append(mohoi)

    stdi = stdstr[ista]
    stdi = float(stdi[0:2]) 
    std.append(stdi)

thick = np.array(thick) 
meanthick = np.mean(thick) 
print('Mean thickness is') 
print(meanthick)

std = np.array(std) 
meanstd = np.mean(std) 
print('Mean std is') 
print(meanstd)

#%% Estimate new thickness, using a formula for expected error from our 2023 paper. 
xiav = 1.08
prc_err = (0.25 * (xiav-1))

newthick = meanthick - meanthick*prc_err
# %%
