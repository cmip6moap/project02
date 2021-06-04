import numpy as np
import os
import sklearn.metrics
from scipy.optimize import curve_fit


def slice_lat(ds):
    return ds.sel(lat=slice(-25, 25))


def ensure_dir(file_path):
    """Check if a directory exists and create it if needed"""
    if not os.path.exists(file_path):
        os.makedirs(file_path)


def days_per_month(month, year):
    """Return the number of days in any month and year"""
    days = [30, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    d = days[int(month)-1]
    if d==28 and int(year)%4==0:
        if int(year)%100==0 and int(year)%400!=0:
            pass
        else:
            d = 29
    return d


def precip_to_mm(ds):
    """Convert precip to mm"""
    if ds.pr.attrs['units']=='kg m-2 s-1':
        ds['pr'] = ds.pr * 24*60**2
        ds.pr.attrs['units']='mm day-1'
    elif ds.pr.attrs['units']=='mm day-1':
        pass
    else:
        raise ValueError('Unrecognised units')
    return ds


def gaus(x,a,x0,sigma):
    """Simple normal distribution function"""
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def fit_gaussian(y, x):
    """Fit a normal gaussian distribution curve to the data. 
    
    Returns 
        [amplitide, mean, width, r^2 statistic]
        4x4 covariance matrix for above values
    """
    popt_f, pcov_f = np.full(4, np.nan, dtype=np.float64), np.full((4,4), np.nan, dtype=np.float64)
    bounds = (np.array([0, -30, 0]), np.array([25, 20, 25]))
    try:
        popt, pcov = curve_fit(gaus,x,y,p0=[8,-5,10], maxfev=8000, bounds=bounds)
        a,x0,sigma = popt
        y_pred = gaus(x,a,x0,sigma)
        r = sklearn.metrics.r2_score(y, y_pred)
        popt_f[:3] = popt
        popt_f[3] = r
        pcov_f[:3, :3] = pcov
    except RuntimeError:
        pass
    return popt_f, pcov_f