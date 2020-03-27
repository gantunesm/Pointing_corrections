import pandas as pd
import numpy as np
import matplotlib as pl

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import get_sun,SkyCoord, EarthLocation, AltAz

# input data
file_data 	= '/home/ccaimapo/scratch/pointing_test/fits.txt'
df_tot 		= pd.read_csv(file_data, sep='\s+\|*\s*', header=None)
N 			= len(df_tot)
act_pos		= (-22.9586, -67.7875, 5190.0)

# define a new data frame with the selected data
df = pd.DataFrame(columns=['Ctime','Freq','Array','Offset X','Error Offset X','Offset Y','Error Offset Y','Modulus Offset','Elevation','Azimuth','Elevation Sun','Azimuth Sun','Angular Distance'])

c = 0
for n in range(N):
	#print(n,'/',N)
	# Get the item 0
	first = df_tot[0][n]
	array = first.split('.')
	array2 = array[2].split(':')
	freq  = int(array2[1].replace('f',''))
	arr_det = array2[0]
	# signal to noise ratio in T, column 13
	snr_T = df_tot[13][n]
	# ctime in days since 2018
	ctime		= float(array[0])
	ctime_days 	= (ctime - 1514764800)/86400

	if arr_det == 'ar4' and freq==150 and snr_T >= 25.0:
		# compute angular distance to the sun
		act_az		= np.radians(float(df_tot[19][n]))
		act_el		= np.radians(float(df_tot[20][n]))
		
		# julian date from ctime
		time		= Time(ctime, format='unix', scale='utc')
		sun_pos		= get_sun(time)
		ACT			= EarthLocation(lat=act_pos[0]*u.deg, lon=act_pos[1]*u.deg, height=act_pos[2]*u.m)
		sun_AzEl	= sun_pos.transform_to(AltAz(obstime=time,location=ACT))
		sun_az		= float(sun_AzEl.az.radian)
		sun_el		= float(sun_AzEl.alt.radian)
		
		#print(sun_az,sun_el,act_az,act_el)
		#exit()
		
		# in radians
		angDist		= 2* np.arcsin(np.sqrt(np.sin(0.5*(sun_el - act_el))**2 + np.cos(sun_el)*np.cos(act_el)*np.sin(0.5*(sun_az - act_az))**2))
		
		deltaX		= float(df_tot[2][n])
		deltaY		= float(df_tot[3][n])
		
		ErrdeltaX		= float(df_tot[4][n])
		ErrdeltaY		= float(df_tot[5][n])

		df.loc[c] = [ctime,freq,arr_det,deltaX,ErrdeltaX,deltaY,ErrdeltaY,np.sqrt(deltaX**2+deltaY**2),act_el,act_az,sun_el,sun_az,angDist]
		c = c + 1

df 	= df.sort_values('Ctime', ascending=True)
df 	= df.drop_duplicates(subset='Ctime', keep='first')

df.reset_index(drop=True,inplace=True)

# save to file
df.to_csv('pointing_corr_arr4_f150_snrT25_FullSeason.csv',index=False)
