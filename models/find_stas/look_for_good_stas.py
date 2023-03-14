#%%
import numpy as np 
import matplotlib.pyplot as plt 
import obspy 
import pandas as pd 
import csv

#%%
fname_sum  = 'ears_sta_summary.csv'
df = pd.read_csv(fname_sum)
# df = df[df.NumEQ >= 200]
# df = df[(df.Long > -120) * (df.Long < -107) * (df.Lat > 25) * (df.Lat < 50)]; 
df = df[(df.Long > -111) * (df.Long < -105) * (df.Lat > 32) * (df.Lat < 40)]; 

# df = df[df['Complexity'] < 0.4]; 
df = df[df['#Net'] == 'XM99']
netkey = '#Net'
print(len(df[netkey].unique()))
# 
# 
fig, ax = plt.subplots(figsize = (20,20))
# handles = []; 
# fact = pd.factorize(df.Network); 


# for iscat = range(length(fact))
#     handles = ax.scatter(df.Long, df.Lat, c=fact[0][iscat], 
#         label = fact[iscat][1] )

handles = []
for i, net in enumerate(df[netkey].unique()):
    handles.append(plt.scatter(df.Long[net==df[netkey]], df.Lat[net==df[netkey]], 
        label = net))

plt.legend()
fig.savefig('temp', dpi = 300)

print(len(df))


# for ista, (net, sta) in enumerate(zip(df['#Net'], df['Station'])):
#     print(net, sta)
#%%
dfsave = netsta = df['#Net'] + ' ' +  df['Station']
# np.savetxt('stationListNew.csv', dfsave.values)
# %%
dfsave.to_csv('stationListNew.csv', index=False, header=False, quoting=csv.QUOTE_NONE,escapechar=",")
# %%
