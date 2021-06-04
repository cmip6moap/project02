import xarray as xr
import numpy as np
import glob
from fit_gaussian_utils import (
    ensure_dir,
    fit_gaussian,
    precip_to_mm,
    slice_lat,
)

root_folder = "/gws/pw/j05/cop26_hackathons/bristol/project02/data/CMIP6histo/zonalAverages"
#output_root_folder = "/gws/pw/j05/cop26_hackathons/bristol/project02/data/CMIP6histo/TRBindices"
output_root_folder = "data/CMIP6histo/TRBindices"


folders = [f for f in glob.glob(root_folder+"/*") if ".sh" not in f]

for folder in folders:
    if 'ymonmean_files' in folder:
        continue
    # create output folder
    output_folder = f"{output_root_folder}/{folder.split('/')[-1]}"
    ensure_dir(output_folder)
    
    # load dataset
    for file in glob.glob(f"{folder}/*.nc"):
        print(f"Processing : {file}")
        output_file = f"{output_folder}/{file.split('/')[-1]}"
        ds = xr.open_dataset(file)
        ds = slice_lat(ds)
        ds = precip_to_mm(ds)

        # run the fitting process
        res = xr.apply_ufunc(fit_gaussian, 
                             ds.pr.astype(np.float64), 
                             ds.lat.astype(np.float64), 
                             input_core_dims=[["lat"], ["lat"]], 
                             output_core_dims=[['popt'], ['popt', 'popt']],
                             vectorize =True)

        # append these to the original data
        ds['gaussian_params'] = res[0]
        ds['gaussian_params_uncertainty'] = res[1]
        ds = ds.assign_coords(popt=['amplitude', 'mean', 'width', 'r2'])

        # select only the fitted parameters
        ds = ds[['gaussian_params', 'gaussian_params_uncertainty']]

        ds.to_netcdf(output_file)
