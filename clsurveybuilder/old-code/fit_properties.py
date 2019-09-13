import numpy as np
from pandas import DataFrame
import sys
import scipy.optimize as opt

def histedges_equalN(x, nbin):
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))



def fit_relation(catalog, nbins=50):
    data = np.zeros((nbins, 4))
    bins = histedges_equalN(catalog['obs_richness'].values, nbins)
    data[:,0] = bins[:-1]
    data[:,1] = bins[1:]
    fitdata = DataFrame(data, columns=['low', 'high', 'median_lamobs', 'median_lnLx'])

    # Bin the data
    for idx, row in fitdata.iterrows():
        subcat = catalog[(catalog['obs_richness'] > row['low']) &
                         (catalog['obs_richness'] <= row['high'])]
        if len(subcat) == 0:
            row['median_lamobs'] = np.NaN
            row['median_lnLx'] = np.NaN
            continue
        row['median_lamobs'] = np.median(subcat['obs_richness'])
        row['median_lnLx'] = np.median(subcat['true_lnLx'])
    fitdata = fitdata.dropna()

    # Set up for fit
    xdata = np.log(fitdata['median_lamobs'].values/30.)
    ydata = fitdata['median_lnLx'].values

    model = lambda x, p1, p2: p1 + p2 * x
    fit, cov = opt.curve_fit(model, xdata, ydata)

    return fitdata, fit

def compute_expectation_bins(cat, xcol, ycol, nbins=30):
    xmin = np.min(cat[xcol])
    xmax = np.max(cat[xcol])
    bins = np.linspace(xmin, xmax, nbins+1)

    data = np.zeros((nbins, 4))
    data[:,0] = bins[:-1]
    data[:,1] = bins[1:]
    fitdata = DataFrame(data, columns=['low', 'high', 'median_lnx', 'median_y'])

    for idx, row in fitdata.iterrows():
        subcat = cat[(cat[xcol] > row['low']) &
                     (cat[xcol] <= row['high'])]
        if len(subcat) == 0:
            row['median_lnx'] = np.NaN
            row['median_y'] = np.NaN
            continue
        row['median_lnx'] = np.median(subcat[xcol])
        row['median_y'] = np.median(subcat[ycol])
    fitdata = fitdata.dropna()

    lnpivot = np.median(cat['lnm200m'])
    xdata = fitdata['median_lnx'] - lnpivot
    ydata = fitdata['median_y']

    model = lambda x, p1, p2: p1*np.exp(x)**p2
    fit, cov = opt.curve_fit(model, xdata, ydata)

    cat['expected_{}'.format(ycol)] = fit[0]*np.exp(cat[xcol] - lnpivot)**fit[1]
    cat['res_{}'.format(ycol)] = cat[ycol] - cat['expected_{}'.format(ycol)]
    return cat, fitdata







