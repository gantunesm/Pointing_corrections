import pandas as pd
import numpy as np
import matplotlib as pl

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import get_sun,SkyCoord, EarthLocation, AltAz

# input data
file_data 	= '/Users/chervias/Documents/Postdoc_FSU/Pointing_correction/fits_s19.txt'
file_name = 'pointing_corr_arr4_f150_snr0_FullSeason19.csv'
df_tot 		= pd.read_csv(file_data, sep='\s+\|*\s*', header=None)
df_tot 		= df_tot.drop_duplicates(subset=0, keep='first')
isSNR = df_tot[15] >= 0.0
df_tot = df_tot[isSNR]
df_tot.reset_index(drop=True,inplace=True)
N 			= len(df_tot)


act_pos		= (-22.9586, -67.7875, 5190.0)
ACT			= EarthLocation(lat=act_pos[0]*u.deg, lon=act_pos[1]*u.deg, height=act_pos[2]*u.m)

# define a new data frame with the selected data
df = pd.DataFrame(columns=['ctime','ctime 2 hr ago','freq','array'])
c = 0
for n in range(N):
	first = df_tot[0][n]
	array = first.split('.')
	array2 = array[2].split(':')
	freq  = int(array2[1].replace('f',''))
	arr_det = array2[0]
	ctime		= float(array[0])
	ctime_2hrago = ctime - 7200.0
	df.loc[c] = [ctime,ctime_2hrago,freq,arr_det]
	c += 1

df.insert(4,'snr',df_tot[15],True)
df.insert(4,'Offset X',df_tot[2],True)
df.insert(4,'Offset Y',df_tot[3],True)
df.insert(4,'Error Offset X',df_tot[4],True)
df.insert(4,'Error Offset Y',df_tot[5],True)
df.insert(4,'Elevation ACT', np.radians(df_tot[20]),True)
df.insert(4,'Azimuth ACT', np.radians(df_tot[19]),True)

# this is the condition for filtering by array
isArray = df['array'] == 'ar4'
# this is the condition for filtering by freq.
isFreq = df['freq'] == 150
df = df[isArray & isFreq ]
df.reset_index(drop=True,inplace=True)

print(df.head())


sun_el_arr = [] ; sun_az_arr = [] ; angdist_arr = [] ; sun_el_2hrago_arr = [] ; sun_az_2hrago_arr = []

for index, row in df.iterrows():
	ctime = row['ctime']
	ctime_2hrago = row['ctime 2 hr ago']
	act_el = row['Elevation ACT']
	act_az = row['Azimuth ACT']
	# julian date from ctime
	time		   = Time(ctime, format='unix', scale='utc')
	time_2hrago	   = Time(ctime_2hrago, format='unix', scale='utc')
	sun_pos		   = get_sun(time)
	sun_pos_2hrago = get_sun(time_2hrago)
	sun_AzEl	   = sun_pos.transform_to(AltAz(obstime=time,location=ACT))
	sun_AzEl_2hago = sun_pos_2hrago.transform_to(AltAz(obstime=time_2hrago,location=ACT))
	sun_az		   = float(sun_AzEl.az.radian)
	sun_el		   = float(sun_AzEl.alt.radian)
	sun_az_2hrago  = float(sun_AzEl_2hago.az.radian)
	sun_el_2hrago  = float(sun_AzEl_2hago.alt.radian)
	angDist		   = 2* np.arcsin(np.sqrt(np.sin(0.5*(sun_el - act_el))**2 + np.cos(sun_el)*np.cos(act_el)*np.sin(0.5*(sun_az - act_az))**2))
	sun_el_arr.append(sun_el);sun_az_arr.append(sun_az) ; angdist_arr.append(angDist) ; sun_el_2hrago_arr.append(sun_el_2hrago) ; sun_az_2hrago_arr.append(sun_az_2hrago)

df.insert(4,'Azimuth Sun',sun_az_arr,True)
df.insert(4,'Elevation Sun',sun_el_arr,True)
df.insert(4,'Angular Dustance',angdist_arr,True)
df.insert(4,'Azimuth Sun 2 hr ago',sun_az_2hrago_arr,True)
df.insert(4,'Elevation Sun 2 hr ago',sun_el_2hrago_arr,True)

# save to file
df.to_csv(file_name,index=False)
