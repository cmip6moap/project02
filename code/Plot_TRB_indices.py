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


filename_TRB_indices='/gws/pw/j05/cop26_hackathons/bristol/project02/data/obs/TRBindices/GPCC/pr_1m_zonAvgNikulin_GPCC_1891-2019_05.nc'


nc_fid = Dataset(filename_TRB_indices,'r')



time = nc_fid.variables['time'][:]
Amplitude = nc_fid.variables['popt'][0]
Mean = nc_fid.variables['popt'][1]
Width = nc_fid.variables['popt'][2]
gauss_params = nc_fid.variables['gaussian_params'][:]



fig = plt.figure(figsize=(10,10))#
ax1 = fig.add_subplot(1,1,1)

plt.plot(time,gauss_params[:,0])  #Amplitude(intensity)
#plt.plot(time,gauss_params[:,1])  #mean
#plt.plot(time,gauss_params[:,0])  #Width


ax1.set_title('TRB indices: Intensities for GPCC')
ax1.set_xlabel('time')




plt.savefig('/gws/pw/j05/cop26_hackathons/bristol/project02/plots/historical/TRB_indices/TRB_index_Amplitude_vrs_time_GPCC.png')
plt.show()
