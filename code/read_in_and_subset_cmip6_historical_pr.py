import os
import numpy as np
import xarray as xr
from datetime import datetime as dt
from pathlib import Path
import iris
import iris.coord_categorisation
from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units


def get_dates(cube, verbose=False):
    """
    Function to get dates from iris cube
    """
    dates = cube.coord('time').units.num2date(cube.coord('time').points)
    dates = [dt(date.year, date.month, date.day) for date in dates]
    if verbose is True:
        print(dates)
    else:
        print(dates[0], '–', dates[-1])
    return(dates)

def make_cmip6_filepath(institute, model, scenario, variant, experiment, table_id, variable, grid, version, time_range,
                        data_root="/badc/cmip6/data/CMIP6/"):
    """
    Make a file path for a cmip6 dataset on JASMIN for a single variable
    """
    # get base path
    path = str(DATA_ROOT / scenario / institute / model / experiment)
    #print(path)
    #print(os.listdir(path))
    
    # get path for variant
    if variant is None:
        # select first variant
        dir_list = os.listdir(path)
        variant_list = [x for x in dir_list if x.startswith('r')]
    else:
        variant_list = [variant]
    
    # update path
    var = [x for x in variant_list if x.startswith('r1i1p1')]
    if len(var) == 0:
        print(variant_list)
        var = [x for x in variant_list if x.startswith('r')]
        path = path + '/' + var[0] + '/' + str(table_id) + '/' + str(variable)
    else:
        path = path + '/' + var[0] + '/' + str(table_id) + '/' + str(variable) 
    #print(path)
    # get path for grid
    if grid is None:
        # select first grid (usually only 1)
        dir_list = os.listdir(path)
        grid_list = [x for x in dir_list if x.startswith('g')]
    else:
        grid_list = [grid]
        
    # update path
    path = path + '/' + str(grid_list[0])
    #print(path)
    
    # get version path
    if version is None:
        dir_list2 = os.listdir(path)
        version_list = [x for x in dir_list2 if x.startswith('v')]
    else:
        version_list = [version]
    
    # update path
    path = path + '/' + str(version_list[0]) + '/'
    print('JASMIN FILEPATH:')
    print(path)
    print('DIRECTORY CONTENTS:')
    print(os.listdir(path))
    return(path+ '*.nc')
    
    

# create dictionary of models and institutes (allows you to loop over models and know the name of the directory that contains the data for that model)
basepath = '/badc/cmip6/data/CMIP6/CMIP/'
institute_list = os.listdir(basepath)
model_inst_dict = {}

# loop over institutes
for inst in institute_list:
    model_list = os.listdir(basepath + inst + '/')
    
    # for each institute list models and store in dictionary
    for model_temp in model_list:
        model_inst_dict[model_temp] = inst
    
    # correction for UKESM which is used by multiple centres - we want MOHC only
    model_inst_dict['UKESM1-0-LL'] = 'MOHC'
    
print(model_inst_dict)



#assert False  # May want to add this if already saved data into dictionaries to prevent re-writing

# Read in precipitation data over domain for CMIP6 models and save data to dictionary

DATA_ROOT = Path("/badc/cmip6/data/CMIP6/")

# dictionary to save subset model output
pr_datasets = {}

# Define region of interest
latmin = -35
latmax = 25
lonmin = 20
lonmax = 30

# variables for JASMIN directory structure
table_id = 'Amon'  # monthly model output
variable_id = 'pr'  # variable code for precipitation in cmip6 model output

# read in monthly data 
# I don't think this is a full list! Need to update
models = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CAMS-CSM1-0', 'CanESM5',
          'CNRM-CM6-1', 'CNRM-ESM2-1', 'FGOALS-f3-L', 'FGOALS-g3', 'HadGEM3-GC31-MM',
          'GISS-E2-1-G', 'INM-CM5-0', 'INM-CM4-8',
          'MPI-ESM1-2-LR', 'NorESM2-LM', 'NorESM2-MM', 'TaiESM1', 'UKESM1-0-LL'] 

# Try for just one model to see if it works
#models = ['HadGEM3-GC31-MM']

# Loop over multiple model experiments and calculate SPEI for all, north and south Ghana
#for expt in ['historical', 'ssp119', 'ssp585']:
for expt in ['historical']:

    for model in models:
        print(model, '', expt.upper())
        institute = model_inst_dict[model]

        if model in ['UKESM1-0-LL']:  #something wrong with UKESM r1i1p1 variant (hdf error)
            variant = 'r2i1p1f2'
        else:
            variant = None

        try:
            # get CMIP6 precip data
            if expt == 'historical':
                scenario = 'CMIP'
            elif expt in ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
                scenario = 'ScenarioMIP'
            
            # get filepath for data for particular model and variable of interest
            fp_hist = make_cmip6_filepath(institute=institute, scenario=scenario, model=model, experiment=expt, variant=variant,
                                          table_id=table_id, variable=variable_id, grid=None, version=None, time_range="*")
            
            # read in data
            pr_data = xr.open_mfdataset(fp_hist)
            pr_data = pr_data.assign_coords(lon=(((pr_data.lon + 180) % 360) - 180)).sortby('lon') # change lons from 0,360 to -180,180
            
            # select data over domain and convert to Iris cube (latmin,latmax, lonmin, lonmax defined above)
            pr_regional_subset = pr_data.sel(lat=slice(latmin,latmax), lon=slice(lonmin,lonmax))
            pr_regional_cube = pr_regional_subset.pr.to_iris()

            #PW: Taking zonal mean            
            pr_regional_cube_zm=pr_regional_cube.collapsed('lon',iris.analysis.MEAN)

            # if you want you can print cube to look at data and print date range
            print()
            print(pr_regional_cube_zm)
            print()
            dates = get_dates(pr_regional_cube_zm)
            pr_datasets[model] = pr_regional_cube_zm


            #PW: setting output file name for each cube
            outpath = '/gws/pw/j05/cop26_hackathons/bristol/project02/data/CMIP6histo/zonalAverages/'+model+'/'
            if not os.path.isdir(outpath):
                os.mkdir(outpath)
            fname = 'pr_1m_zonAvgNikulin_'+model+'_'+expt+'_s'+str(dates[0].year)+'.nc'
            print('SAVING TO:', outpath + fname)
            iris.save(pr_regional_cube_zm, outpath + fname)
            print()

        except FileNotFoundError:
            print(model, ' has no ' + expt.upper() + ' output')
            continue
        
    """# change output directory to somewhere you can save
    outpath = '/home/users/jcabaker/bristol_cmip6_hack/save_files/'
    fname = 'cmip6_' + expt + '_pr_dict.npy'
    print('SAVING TO:', outpath + fname)
    print()
    #np.save(outpath + fname, pr_datasets)  # uncomment to save"""
