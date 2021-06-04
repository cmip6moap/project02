#Michael Baidu CMIP6 Hackathon June 2021

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp1
from netCDF4 import Dataset
#------------------------------------------
import iris
import iris.plot as iplt
import iris.coords
from iris.coords import DimCoord
from iris.coords import Coord
#-----------------------------------------------
import os

# New package for this week
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#-----------------------------------------------
import cartopy
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Set the font dictionaries (for plot title and axis titles)
title_font = {'fontname':'Arial', 'size':'11', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more spac
quiverkey_font = {'color':'black'} # Bottom vertical alignment for more spac
plt.rcParams.update({'font.size': 11})


models = ['GPCC','ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CAMS-CSM1-0', 'CanESM5','CNRM-CM6-1','CNRM-ESM2-1','FGOALS-g3','HadGEM3-GC31-MM','GISS-E2-1-G','INM-CM5-0', 'INM-CM4-8','MPI-ESM1-2-LR', 'NorESM2-LM', 'NorESM2-MM', 'TaiESM1', 'UKESM1-0-LL']

n=0
fig = plt.figure(figsize=(10,10))#
for model in models:
	n=n+1
	if model =='GPCC':
		data_type='obs'
		period='1891-2019_05'
	else:
		data_type='CMIP6histo'
		period='historical_s1850'
	

	filename_TRB_indices='/gws/pw/j05/cop26_hackathons/bristol/project02/data/'+str(data_type)+'/TRBindices/'+str(model)+'/pr_1m_zonAvgNikulin_'+str(model)+'_'+str(period)+'.nc'

	

	nc_fid = Dataset(filename_TRB_indices,'r')



	time = nc_fid.variables['time'][:]
	Amplitude = nc_fid.variables['popt'][0]
	Mean = nc_fid.variables['popt'][1]
	Width = nc_fid.variables['popt'][2]
	gauss_params = nc_fid.variables['gaussian_params'][:]
	data=gauss_params
	
	#seas_cyc=np.array([12])
	#for i in range(12):
		#seas_cyc[i]=np.nanmean(data[i::12,0])



	
	
	ax1 = fig.add_subplot(3,6,n)

	plt.plot(time,gauss_params[:,0])  #Amplitude(intensity)
	#plt.plot(time,seas_cyc[:,0])  #Amplitude(intensity)
	#plt.plot(time,gauss_params[:,1])  #mean
	#plt.plot(time,gauss_params[:,0])  #Width


	ax1.set_title('Intensities- '+str(model)+'')
	ax1.set_xlabel('time')




plt.savefig('/gws/pw/j05/cop26_hackathons/bristol/project02/plots/historical/TRB_indices/TRB_index_Amplitude_vrs_time_'+str(model)+'_'+str(period)+'.png')
plt.show()

	
	
	
	
	
	
	
	
	
	
	
	
