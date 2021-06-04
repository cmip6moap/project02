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
title_font = {'fontname':'Arial', 'size':'12', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more spac
quiverkey_font = {'color':'black'} # Bottom vertical alignment for more spac
plt.rcParams.update({'font.size': 12})

models = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CAMS-CSM1-0', 'CanESM5','CNRM-CM6-1','CNRM-ESM2-1','FGOALS-g3','HadGEM3-GC31-MM','GISS-E2-1-G','INM-CM5-0', 'INM-CM4-8','MPI-ESM1-2-LR', 'NorESM2-LM', 'NorESM2-MM', 'TaiESM1', 'UKESM1-0-LL']


for model in models:

	filename_CMIP6_zonal_mean='/gws/pw/j05/cop26_hackathons/bristol/project02/data/CMIP6histo/zonalAverages/ymonmean_files/pr_1m_zonAvgNikulin_'+str(model)+'_historical_s1850_ymonmean.nc'




	nc_fid = Dataset(filename_CMIP6_zonal_mean,'r') 

	slats = nc_fid.variables['lat'][:]
	#slons = nc_fid.variables['lon'][:]
	time = nc_fid.variables['time'][:]

	plot_lats = slats

	precip=nc_fid.variables['pr'][:]

	#precip(time, lat) 

	fig = plt.figure(figsize=(10,10))#
	#ax1 = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
	ax1 = fig.add_subplot(1,1,1)

	nx=np.arange(np.size(time))
	print(nx)


	precip = np.transpose(precip)



	im1=ax1.contourf(nx,plot_lats[:],precip[:,:],np.arange(0.00002,0.00012+0.00002,0.00002),cmap=plt.cm.bwr,extend='both')
	ax1.set_xticks(nx)
	ax1.set_xticklabels(('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))
	ax1.set_ylabel('latitude $^{o}$')
	cax = fig.add_axes([0.1,0.055,0.8,0.02]) #left, bottom, width, height
	cbar = fig.colorbar(im1, cax=cax, extend='both',orientation='horizontal')
	ax1.set_title('Zonal mean '+str(model)+' historical precip (1981-2014)')
	cbar.ax.set_xlabel('zonal mean precip (mm/month)')


	#ax1.clabel(im6, inline=1,fmt='%1.2f', fontsize=10)
	plt.savefig('/gws/pw/j05/cop26_hackathons/bristol/project02/plots/historical/Zonal_mean_'+str(model)+'_historical_precip.png')
	plt.show()
