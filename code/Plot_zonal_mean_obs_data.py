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
title_font = {'fontname':'Arial', 'size':'16', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more spac
quiverkey_font = {'color':'black'} # Bottom vertical alignment for more spac
plt.rcParams.update({'font.size': 18})


#filename_GPCC_zonal_mean='/gws/pw/j05/cop26_hackathons/bristol/project02/data/obs/zonalAverages/GPCC/pr_1m_zonAvgNikulin_GPCC_1891-2019_05_ymonmean.nc'
filename_GPCC_zonal_mean='/gws/pw/j05/cop26_hackathons/bristol/project02/data/obs/zonalAverages/pr_1m_zonAvgNikulin_GPCC_1981-2014_05_ymonmean.nc'



nc_fid = Dataset(filename_GPCC_zonal_mean,'r') 

slats = nc_fid.variables['lat'][:]
#slons = nc_fid.variables['lon'][:]
time = nc_fid.variables['time'][:]

plot_lats = slats

precip=nc_fid.variables['precip'][:]

#precip(time, lat) 

fig = plt.figure(figsize=(10,10))#
#ax1 = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
ax1 = fig.add_subplot(1,1,1)

nx=np.arange(np.size(time))
print(nx)


precip = np.transpose(precip)



im1=ax1.contourf(nx,plot_lats[:],precip[:,:],cmap=plt.cm.bwr,extend='both')
ax1.set_xticks(nx)
ax1.set_xticklabels(('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))
cax = fig.add_axes([0.1,0.04,0.8,0.02]) #left, bottom, width, height
cbar = fig.colorbar(im1, cax=cax, extend='both',orientation='horizontal')
ax1.set_title('Zonal mean GPCC precip (1981-2019)')
cbar.ax.set_xlabel('zonal mean precip (mm/month)')


#ax1.clabel(im6, inline=1,fmt='%1.2f', fontsize=10)
plt.savefig('/gws/pw/j05/cop26_hackathons/bristol/project02/plots/historical/Zonal_mean_GPCC_precip.png')
plt.show()
