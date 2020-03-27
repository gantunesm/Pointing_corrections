#!/usr/bin/env python
# coding: utf-8

# In[5]:


import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import seaborn as sns
from scipy.stats import pearsonr
from moby2.aux_data import hk
import moby2
import zipfile
import os, sys, glob
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import get_sun,SkyCoord, EarthLocation, AltAz


# In[3]:


# only this part needs to be modified
# Change the obs number and the ctime limits
obs = 0
# Change the data base created before
df = pd.read_csv('pointing_corr_arr4_f150_snrT25_FullSeason.csv')
# this is the name of the output data base
output_db = 'pointing_corr_arr4_f150_snrT25_obs%i_FullSeason_withEverything.csv'%obs

# obs 0
t0 = 1525652963.0
t1 = 1525962718.0

# obs 1
#t0 = 1527281490.
#t1 = 1528231280.

# obs 2
#t0 = 1536609972.
#t1 = 1544458086.


# In[6]:


getter_apex = moby2.aux_data.apex.Temperature()
names = ['Enc_El_Deg', 'az_mean']
getters_boresight = [hk.HKChannel(n) for n in names]
act_pos		= (-22.9586, -67.7875, 5190.0)
ACT			= EarthLocation(lat=act_pos[0]*u.deg, lon=act_pos[1]*u.deg, height=act_pos[2]*u.m)
N 			= len(df)

# Find all relevant hk files.
# the prefixes go between 
prefixes = np.arange(int(str(t0)[0:5]),int(str(t1)[0:5])+1,1)


# In[7]:


roots = [x['params']['root']
         for x in moby2.user_cfg['filebase']['components']]
hk_files = {}
for p in prefixes:
    _hkf = []
    for r in roots:
        _hkf.extend(glob.glob('%s/%s/*.hk.zip' % (r, p)))
    for x in _hkf:
        hk_files[os.path.split(x)[1]] = x
# Field study.
keys = sorted(hk_files.keys())

#FM = moby2.DirfileManager(hk_files[keys[1]])

pri_fields = []

for key in keys:
	#print(key)
	FM = moby2.DirfileManager(hk_files[key])
	try:
		for line in FM.get_fields():
			w = line[1:]
			if len(w) == 0 or w[0][0] == '#':
				continue
			if w[0].startswith('T_Pan'):
				if w[0] not in pri_fields:
					pri_fields.append(w[0])
	except UnicodeDecodeError:
		continue

#pri_fields.remove('T_Pan_R1_2')
FIELD_SETS = [
    ('pri', pri_fields, 4)]

print(FIELD_SETS)

# Load all that data, keeping times we care about.
data = {'t': []}
for k in keys:
    #print(k)
    FM = moby2.DirfileManager(hk_files[k])
    i = FM.get_frame_count('cpu_s')
    t = FM.load_channel('cpu_s', 0, i.n_samples)[::400].astype('float')
    s = (t >= t0) * (t < t1)
    if not np.any(s): continue
    #az  = FM.load_channel('Enc_Az_Deg_Astro', 0, i.n_samples)[::400].astype('float')
    #alt = FM.load_channel('Enc_El_Deg', 0, i.n_samples)[::400].astype('float')
    data['t'].append(t[s])
    #print(t[s])
    #data['az'].append(az[s])
    #data['alt'].append(alt[s])
    for (mir, fields, ref_i) in FIELD_SETS:
        for f in fields:
            i = FM.get_frame_count(f)
            try: 
                z = FM.load_channel(f, 0, i.n_samples)
            except RuntimeError:
                # skip that field
                print(data['t'][-1][0])
                exit()
            z = z[::i.spf]
            #print f, i.n_samples, z.shape, i.spf
            if not f in data:
                data[f] = []
            data[f].append(z[s])
        
for k,v in list(data.items()):
    data[k] = np.hstack(v)


# In[ ]:


keys = sorted(data.keys())
nkeys = len(keys)

# this will select from the full db between the ctimes of interest
condition = np.logical_and(t0 <= df['Ctime'],df['Ctime'] <= t1)

df = df[condition]
df.reset_index(inplace=True,drop=True)
N 			= len(df)
temp = np.zeros((N,nkeys))

c = 0
temp_array = []
temp_2hrago_arr = []
el_2hrago_arr = []
az_2hrago_arr = []
sun_el_2hrago_arr = []
sun_az_2hrago_arr = []
angdistance_2hrago_arr = []

for n in range(N):
	ctime        = df['Ctime'][n]
	ctime_2hrago = ctime - 7200.0
	# we assume 1 TOD is 11 minutes long
	ctime_f = ctime + 660.0
	ctime_f_2hrago = ctime_2hrago + 660.0
	#print(ctime,ctime_f)
	nsamples, Temp, std_temp = getter_apex.get_average(ctime,ctime_f)
	nsamples, Temp_2hrago, std_temp = getter_apex.get_average(ctime_2hrago,ctime_f_2hrago)
	temp_array.append(Temp)
	temp_2hrago_arr.append(Temp_2hrago)
	i_ctime = (np.abs(data['t']-ctime)).argmin()
	i_ctime_f = (np.abs(data['t']-ctime_f)).argmin()
	for nn in range(nkeys):
		temp[n,nn] = np.mean(data[keys[nn]][i_ctime:i_ctime_f])
	# look for values 2 hours ago
	boresight_2hrago = [g.get_nearest(ctime_2hrago) for g in getters_boresight]
	el_2hrago    = np.radians(boresight_2hrago[0][0])
	az_2hrago    = np.radians(boresight_2hrago[1][0])
	el_2hrago_arr.append(el_2hrago)
	az_2hrago_arr.append(az_2hrago)
	# julian date from ctime
	time_2hrago	  = Time(ctime_2hrago, format='unix', scale='utc')
	sun_pos_2hrago= get_sun(time_2hrago)
	sun_AzEl	  = sun_pos_2hrago.transform_to(AltAz(obstime=time_2hrago,location=ACT))
	sun_az_2hrago = float(sun_AzEl.az.radian)
	sun_el_2hrago = float(sun_AzEl.alt.radian)
	sun_el_2hrago_arr.append(sun_el_2hrago)
	sun_az_2hrago_arr.append(sun_az_2hrago)
	angDist_2hrago= 2*np.arcsin(np.sqrt(np.sin(0.5*(sun_el_2hrago - el_2hrago))**2 + np.cos(sun_el_2hrago)*np.cos(el_2hrago)*np.sin(0.5*(sun_az_2hrago - az_2hrago))**2))
	angdistance_2hrago_arr.append(angDist_2hrago)

# insert columns in the data frame		
for nn in range(nkeys-1):
	df.insert(10,keys[nn],temp[:,nn],True)	
df.insert(6,"APEX T", temp_array, True)
df.insert(7,"APEX T 2 hr ago",temp_2hrago_arr,True)
df.insert(8,"Elevation 2 hr ago",el_2hrago_arr,True)
df.insert(9,"Azimuth 2 hr ago",az_2hrago_arr,True)
df.insert(10,"Angular Distance 2 hr ago",angdistance_2hrago_arr,True)
df.insert(11,"Azimuth Sun 2 hr ago",sun_az_2hrago_arr,True)
df.insert(12,"Elevation Sun 2 hr ago",sun_el_2hrago_arr,True)

df.to_csv(output_db,index=False)


# In[ ]:




