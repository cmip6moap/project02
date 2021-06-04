#Script to plot zonal-mean rainfall against latitude and month of the year, averaged over 1981-2014

import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np

import iris
import iris.coord_categorisation  
import math
import xarray as xr

#Python function to return the index of a numpy array where the array value is closest to "value".
def closest(array,value):
  ind=abs(array-value).argmin()
  return ind

#James' code to use with converting GPCC values to mm/day below
days = [30, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
def days_per_month(month, year):
    d = days[int(month)-1]
    if d==28 and int(year)%4==0:
        if int(year)%100==0 and int(year)%400!=0:
            pass
        else:
            d = 29
    return d


#Function to return seasonal cycle computed from a time series
def get_iris_seascyc(cube):

    cube2=cube.copy()
    coord_names=[coord.name() for coord in cube2.coords()]

    #Get the average over all years
    if 'month' not in coord_names:
        iris.coord_categorisation.add_month(cube2,'time',name='month')
    cube_sc=cube2.aggregated_by('month', iris.analysis.MEAN)
    
    return cube_sc


#Function to return levels for a contour plot, with nlevels between the max and min data values. More levels may be returned if shift and centred_shift options are used, to make sure shifted levels cover the data range.
#By default, if the data includes +ve and -ve values, the levels are shifted so one of them is zero. Set shift to a different value to ensure that the levels are shifted so that one of them equals that value. Set the value of shift to fall half way between two levels with centred_shift option.
#Set equal_pos_neg_levs to have the positive and negative levels be mirrors of each other.
#Set log to get logarithmic levels (base 10).
#Set log_round_min_max to, when log=True, round the min and max level values to integers.
#Set clip_outliers to restrict the maximum and minimum levels to lie no further than clip_sd standard deviations from the median, where the standard deviation is calculated from the data clipped so that the values lie in the range of the middle 95% of the original.
def get_levels(data,nlevels=9, shift=None, centred_shift=False, equal_pos_neg_levs=False, log=False, log_round_min_max=True, clip_outliers=False, clip_sd=3):

    if type(data)==list:
        data=[x for x in data if x is not None]  #excluding None elements

    if type(data) not in [np.ndarray, np.ma.core.MaskedArray]:
        try:
            data=np.ma.concatenate(data) #it seems that this needs to be done rather than np.array() if data is a list of masked arrays. I think it should yield the same results in other cases.
        except:
            data=np.array(data)
    data=data[(-np.inf<data) & (data<np.inf) & (data is not np.nan)]  #ignore non-finite values
    assert np.min(data)<np.max(data)  #check that not all values are equal

    if log:
        data[data<=0]=np.min(data[data>0])  #set minimum of data to minimum positive value
        data=np.log10(data)

    if clip_outliers:
        #Get range of levels ignoring values more than clip_sd s.d. from the median (with data being clipped to lie in the range of the middle 95% of the original) if specified
        median=np.ma.median(data)  #use np.ma.median so this works with masked data 
        data_clipped=copy.deepcopy(data)
        data_clipped[data<np.percentile(data,2.5)]=np.percentile(data,2.5)  #using percentiles rather than indexing a sorted list of the array values allows this to work even when data is a small array. (I've not checked if this works for masked data.)
        data_clipped[data>np.percentile(data,97.5)]=np.percentile(data,97.5)
        if type(data)==np.ma.core.MaskedArray:  #for masked arrays, undo any changes to the mask from above
            data_clipped.mask=data.mask
        std=np.ma.std(data_clipped)
        levmax=np.nanmin([ np.max(data), abs(median)+clip_sd*std ])
        levmin=np.nanmax([ np.min(data), -abs(median)-clip_sd*std ])
        
    else:
        levmin=round_to_n(np.min(data),2,round_down=1)
        levmax=np.max(data)

    
    if shift and shift<levmin:
        levmin=shift
    
    if shift and shift>levmax:
        levmax=shift
        
    #For log levels, if log_round_min_max=True, make sure the min and max are integers. If the resulting step between levels is less than one, then choose the step so the number of levels is as close to nlevels as possible whilst the levels include the integer values between levmin and levmax.
    if log and log_round_min_max:
        levmin=np.floor(levmin)
        levmax=np.ceil(levmax)
        step0=(levmax - levmin) / (nlevels-1)
        if step0<1:
            step=1./round(1./step0)
            nlevels=1+(levmax - levmin) / step
        else:
            step=step0
        levels = round_to_n(np.arange(nlevels) * step + levmin, 2)  #round levels rather than step so that integer values stay as integers.
    else:
        step = round_to_n((levmax - levmin) / (nlevels-1), 2, round_up=1)
        levels = np.arange(nlevels) * step + levmin
    #print levmin, levmax, step, levels

        if shift is None and levmax>0 and levmin<0:
            shift=0
    
    #if levmax>shift and levmin<shift: 
    if shift is not None:
        if centred_shift:  #set the value of shift to fall halfway between two levels
            first_lev_over=levels[levels>=shift][0]
            first_lev_under=levels[levels<=shift][-1]
            #print first_lev_over,first_lev_under,((first_lev_over+first_lev_under)/2.-shift)
            levels=levels-((first_lev_over+first_lev_under)/2.-shift)
        else:  #shift levels so one of them equals shift
            shift_ind=closest(levels,shift)
            #print shift_ind, levels[shift_ind], levels[shift_ind]-shift
            levels=levels-(levels[shift_ind]-shift)
        
        #Add an extra level above or below if the shifted levels do not cover the data range any more
        if levmin<np.min(levels):
            levels=np.concatenate(([np.min(levels)-step],levels))
        if levmax>np.max(levels):
            levels=np.concatenate((levels,[np.max(levels)+step]))
    
    if equal_pos_neg_levs:
        npos=len(levels[levels>0])
        nneg=len(levels[levels<0])
        if npos > nneg:
            levels=np.concatenate((-levels[levels>=0][::-1], levels[levels>=0]))
        elif nneg > npos:
            levels=np.concatenate((levels[levels<=0], -levels[levels<=0][::-1]))
        levels=np.unique(levels)  #remove duplicate zeroes if applicable

    return levels


#Class for getting colourbar to be centred at a given value. From https://matplotlib.org/3.2.2/gallery/userdemo/colormap_normalizations_custom.html.
class MidpointNormalize(matplotlib.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

    
#Function to round numbers to a given number n of significant figures
#Modified from Roy Hyunjin Han's answer at http://stackoverflow.com/questions/3410976/how-to-round-a-number-to-significant-figures-in-python
def round_to_n(x,n,round_up=None,round_down=None):

    from numbers import Number

    if round_up and round_down:
        print('round_to_n: Cannot set both round_up and round_down')
        sys.exit()

    if isinstance(x,Number):  #for number arguments, convert to an array, and convert back at the end.
        num_flag=True
        x=np.array([x])
    else:
        num_flag=False
    
    if type(x)==list:
        x=np.array(x)

    x_shape=x.shape
    x=np.ravel(x)
    ret=np.zeros(len(x))
    for ind,i in enumerate(x):
        if not np.isfinite(i):  #to allow infinite or NaN values to be handled
            ret[ind]=i
        elif i==0:
            ret[ind]=0
        else:
            tens=int(np.floor(np.log10(abs(i)))) - (n - 1)
            ret[ind]=round(i, -tens)

            if round_up and ret[ind]<i:
                ret[ind]+=math.pow(10,tens)
            if round_down and ret[ind]>i:
                ret[ind]-=math.pow(10,tens)
    
    ret=ret.reshape(x_shape)

    if num_flag:
        ret=ret[0]
    
    return ret





datasets = ['GPCC', 'CHIRPS', 'MSWEP', 'ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CAMS-CSM1-0', 'CanESM5',
          'CNRM-CM6-1', 'CNRM-ESM2-1', 'FGOALS-f3-L', 'FGOALS-g3', 'HadGEM3-GC31-MM',
          'GISS-E2-1-G', 'INM-CM5-0', 'INM-CM4-8',
          'MPI-ESM1-2-LR', 'NorESM2-LM', 'NorESM2-MM', 'TaiESM1', 'UKESM1-0-LL'] 

#datasets=['GPCC','HadGEM3-GC31-MM']  #testing

year_start=1981
year_end=2014

obs_datasets=['GPCC', 'CHIRPS', 'MSWEP']  #Used to set file paths correctly

#Loading data
cube_dict={}
for d in datasets:

    path='/gws/pw/j05/cop26_hackathons/bristol/project02/data/'
    if d in obs_datasets:
        path+='obs/'
    else:
        path+='CMIP6histo/'
    path+='zonalAverages/'+d+'/'
    
    f='pr_1m_zonAvgNikulin_'
    if d=='GPCC':
        f+='GPCC_1891-2019_05'
        cube_name='precip'
    else:
        f+='*'
        cube_name=None
    f+='.nc'
    
    cube=iris.load_cube(path+f, cube_name)
    
    t_constraint=lambda x: year_start<=x.point.year<=year_end
    constraint=iris.Constraint(coord_values={'time':t_constraint})
    cube=cube.extract(constraint)

    #Get cube units in mm/day
    unit=cube.units.symbol
    if d in obs_datasets:  #convert from monthly totals to mm/day, by converting to xarray DataArray, converting units based on James' code and converting back to an iris cube
        da=xr.DataArray.from_iris(cube)
        days = xr.apply_ufunc(days_per_month, da.time.dt.month, da.time.dt.year, vectorize=True)
        #ds['pr'] = ds.pr / days
        pr = da.values / np.tile(days.values[:,np.newaxis], (1,da.values.shape[1]))
        da.values = pr
        #da.attrs["units"] = 'mm/day'
        cube=da.to_iris()
    elif np.any([x in unit for x in ['kg m-2 s-1','kg m**-2 s-1','m-2.kg.s-1','mm s-1']]):
        cube.data *= 3600*24
    elif unit=='m.s-1':
        cube2.data *= 1000*3600*24
    elif unit=='2.77777777777778e-07 m.s-1':  #mm/hr
        cube2.data *= 24
    
    cube.units='mm/day'

    cube_dict[d]=cube
    


#Get climatological seasonal cycle, and ratio with respect to the first dataset. This requires regridding the data - use the highest resolution grid, so that linear interpolation will be conservative of the total quantity (I think). Identify this grid by finding the latitude axis with the highest number of points between 25S-25N.
cube_lats_list=[cube_dict[d].coord('latitude').points for d in datasets]
cube_lats_list_restricted=[lats[np.where((-25<=lats) & (lats<=25))] for lats in cube_lats_list]
max_lats_no=np.max([len(lats) for lats in cube_lats_list_restricted])
lats_interp=[lats for lats in cube_lats_list_restricted if len(lats)==max_lats_no][0]

cube_sc_dict={}
cube_sc_regrid_dict={}
sc_ratio_dict={}

#Get seasonal cycle for first dataset and regrid it
cube_sc0=get_iris_seascyc(cube_dict[datasets[0]])
cube_sc_dict[datasets[0]]=cube_sc0
cube_sc_regrid0=cube_sc0.interpolate([('latitude',lats_interp)], iris.analysis.Linear())
cube_sc_regrid_dict[datasets[0]]=cube_sc_regrid0

for ind,d in enumerate(datasets):
    cube_sc_dict[d]=get_iris_seascyc(cube_dict[d])

    if ind>0:
        #Calculate ratio with respect to first dataset
        
        #Regrid cube to share the grid of the first dataset
        #Using linear interpolation - I'm not completely sure this will conserve total rainfall - perhaps not if any model has a much finer resolution than the first dataset
        cube_sc=cube_sc_dict[d]
        cube_sc_regrid=cube_sc.interpolate([('latitude',lats_interp)], iris.analysis.Linear())
        
        #Mask where rainfall in both this and the first dataset is less than 1mm/day, to avoid showing points where the ratio is very large because the denominator is very small without there being a large error
        cube_sc_regrid.data.mask[(cube_sc_regrid.data < 1) & (cube_sc_regrid_dict[datasets[0]].data < 1)]=True
        
        cube_sc_regrid_dict[d]=cube_sc_regrid
        sc_ratio_dict[d]=cube_sc_regrid.data/cube_sc_regrid_dict[datasets[0]].data


#PLOTTING###################################################################################################
        
#Plotting the seasonal cycle in each dataset
if len(datasets)==1:
    nrows=1
elif 2<=len(datasets)<=6:
    nrows=2
elif 7<=len(datasets)<=9 or len(datasets)==21:
    nrows=3
else:
    nrows=4

ncolumns=(len(datasets)-1)//nrows+1

fig,axarr=plt.subplots(nrows=nrows,ncols=ncolumns, figsize=(12.6,6))
if ncolumns==1:  
    axarr=axarr[:,np.newaxis]
if nrows==1:  
    axarr=axarr[np.newaxis,:]
    
levels=get_levels([cube_sc_dict[d].data.ravel() for d in datasets])

month_axis=range(1,13)

for ind,d in enumerate(datasets):
    row_ind=ind//ncolumns
    col_ind=ind % ncolumns
    plt.sca(axarr[row_ind, col_ind])
    lat_axis=cube_sc_dict[d].coord('latitude').points
    plt.contourf(month_axis, lat_axis, cube_sc_dict[d].data.T, vmin=levels.min(), vmax=levels.max(), levels=levels, extend='both')
    
    #testing that interpolated data looks reasonable
    """lat_axis=lats_interp
    plt.contourf(month_axis, lat_axis, cube_sc_regrid_dict[d].data.T, vmin=levels.min(), vmax=levels.max(), levels=levels, extend='both')"""
    
    plt.title(d)
    
    #Put axis labels on edge plots only
    if row_ind==nrows-1:
        plt.xlabel('Month')
    if col_ind==0:
        plt.ylabel('Latitude')

    #Put colour bar on rightmost plots only
    if col_ind==ncolumns-1:
        plt.colorbar()

    plt.ylim(-25,25)
        
plt.tight_layout()
#Command for saving - probably best to comment it here and run it manually if it's wanted.
#plt.savefig('/gws/pw/j05/cop26_hackathons/bristol/project02/plots/historical/plot_cmip6_historical_pr_zm_seascyc.png')


#Plotting ratios of values in each dataset to those of the first dataset, to indicate the size of biases
if len(datasets)-1==1:
    nrows=1
elif 2<=len(datasets)-1<=6:
    nrows=2
elif 7<=len(datasets)-1<=9 or len(datasets)-1==21:
    nrows=3
else:
    nrows=4

ncolumns=(len(datasets)-1-1)//nrows+1

fig_ratio,axarr_ratio=plt.subplots(nrows=nrows,ncols=ncolumns, figsize=(12.6,6))
if ncolumns==1:  
    axarr_ratio=axarr_ratio[:,np.newaxis]
if nrows==1:  
    axarr_ratio=axarr_ratio[np.newaxis,:]
    
levels_ratio=np.arange(0.5,1.55,0.1)
lat_axis=lats_interp

cmap='BrBG'
norm=MidpointNormalize(midpoint=1,vmin=levels_ratio.min(),vmax=levels_ratio.max())

for ind,d in enumerate(datasets[1:]):
    row_ind=ind//ncolumns
    col_ind=ind % ncolumns
    plt.sca(axarr_ratio[row_ind, col_ind])
    plt.contourf(month_axis, lat_axis, sc_ratio_dict[d].T, vmin=levels_ratio.min(), vmax=levels_ratio.max(), levels=levels_ratio, extend='both', cmap=cmap)
    
    #Put colour bar on rightmost plots only
    if col_ind==ncolumns-1:
        plt.colorbar()

    plt.title(d)
    
    #Put axis labels on edge plots only
    if row_ind==nrows-1:
        plt.xlabel('Month')
    if col_ind==0:
        plt.ylabel('Latitude')

    plt.ylim(-25,25)
        
plt.tight_layout()

#Command for saving - probably best to comment it here and run it manually if it's wanted.
#plt.savefig('/gws/pw/j05/cop26_hackathons/bristol/project02/plots/historical/plot_cmip6_historical_pr_zm_seascyc_bias.png')

plt.show()

