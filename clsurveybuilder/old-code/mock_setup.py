import numpy as np
import hmf_library as hmflib
from scipy.interpolate import interp1d, interp2d
import os, sys
from scipy.integrate import quad


def build_splines(cfgin, sinfo):
    '''Build the splines for the mass conversion once'''

    # Build the mass conversion spline
    lnm200m_full, _ = hmflib.load_tinker_hmf(0.1, sinfo, -1., 17.)
    lnm200m = lnm200m_full[::4]

    zmin = cfgin['SurveyInfo'].getfloat('zlo')
    zmax = cfgin['SurveyInfo'].getfloat('zhi')
    zgrid = np.linspace(zmin, zmax, 10)
    lnm500c = np.zeros((len(zgrid), len(lnm200m)))

    for i in range(len(zgrid)):
        lnm500c[i,:] = hmflib.mass_conversion(lnm200m, sinfo, zgrid[i])
    mass_conversion = interp2d(lnm200m, zgrid, lnm500c, kind='cubic', bounds_error=True)

    # P(r) \propto 1 - r^2 spline
    rprior = lambda r: 1.0 - r*r
    rlist = np.linspace(-1.0, 1.0, 501)

    newparts = np.vectorize(quad)(rprior, rlist[:-1], rlist[1:], limit=50)[0]
    cdf = np.zeros(1+len(newparts))
    for i in range(1, len(cdf)):
        cdf[i] = cdf[i-1] + newparts[i-1]
    cdf = cdf/cdf[-1]
    pr_spline = interp1d(cdf, rlist, kind='cubic')


    # Return the both in a dictionary
    return {'mass_conversion' : mass_conversion, 'P(r)' : pr_spline}


# =================================================================================================
def build_file_names(cfgin, run_number, modifier=''):
    '''Make a dictionary with all of the required output file names.

    Parameters
    ----------
    cfgin : ConfigParser(object)
        Inputs to the mcmc from the ini config file
    run_number : int
        Index for mock realization
    modifier : str
        String to append to the end of each file name

    Returns
    -------
    file_dict : dict
        Dictionary of each output file name
    '''
    build_directories(cfgin)

    try:
        root = os.path.join(os.environ['DATADIR'], "xray-scatter/mocks/{}/".format(cfgin['General']['synth_catalog']))
    except KeyError:
        root = "mocks/{}/".format(cfgin['General']['synth_catalog'])

    if modifier == '':
        modifier = str(run_number)

    file_dict = {"model" : "{}mock_model_{}.ini".format(root, modifier),
                 "mass_catalog" : "{}mass_catalogs/synth_cat_{}".format(root, modifier),
                 "obs_catalog" : "{}obs_catalogs/prop_cat_{}".format(root, modifier),
                 "richest" : "{}richest{}_{}.dat".format(root, cfgin['SurveyInfo'].getint('nclusters'), modifier)}

    return file_dict


def import_cosmology(cfgin):
    sinfo = hmflib.SurveyInfo(cfgin)
    cosmology = {'flat': True, 'H0': sinfo.H0, 'Om0': sinfo.OmM, 'Ob0': sinfo.Omb,
                 'sigma8': sinfo.sigma8, 'ns': sinfo.ns}
    return sinfo, cosmology


def build_model_dict(model):
    """Make a dictionary of model parameters

    Parameters
    ----------
    model : ConfigParser(object)
        Input fiducial model

    Returns
    -------
    params : dict
        Fiducial model by parameter name
    """
    params = {'r' : model['Model'].getfloat('r'),
              'mgas0' : model['Model'].getfloat('mgas0'),
              'alphamg' : model['Model'].getfloat('alphamg'),
              'sigma0mg' : model['Model'].getfloat('s0mg'),
              'lam0' : model['Model'].getfloat('lam0'),
              'alphalam' : model['Model'].getfloat('alphalam'),
              'sigma0lam': model['Model'].getfloat('s0lam'),
              'shmf': model['Model'].getfloat('shmf'),
              'qhmf': model['Model'].getfloat('qhmf')}

    return params


def save_observed_subcatalog(catalog, filedict):
    subcat = catalog[['lnm200m', 'true_richness', 'obs_richness', 'true_lnLx']]
    subcat.to_csv(filedict['obs_catalog'], sep='\t', index_label='id', float_format='%.6e')
    return


def save_observed_catalog(catalog, filedict):
    """Save the cluster catalog to a file

    Parameters
    ----------
    catalog : DataFrame
        The full cluster catalog
    fildict : dict
        Dictionary of file names
    """
    subcat = catalog[['obs_richness', 'obs_richness_std', 'obs_mgas', 'obs_mgas_std', 'lnm200m']]
    subcat.to_csv(filedict['obs_catalog'], sep='\t', index_label='id', float_format='%.6e')


def build_directories(cfgin):
    try:
        root = os.path.join(os.environ['DATADIR'], "xray-scatter/mocks/".format(cfgin['General']['synth_catalog']))
    except KeyError:
        root = "mocks/".format(cfgin['General']['synth_catalog'])

    pathlist = [os.path.join(root, cfgin['General']['synth_catalog']),
                os.path.join(root, cfgin['General']['synth_catalog'], 'mass_catalogs'),
                os.path.join(root, cfgin['General']['synth_catalog'], 'obs_catalogs')]

    for path in pathlist:
        if not os.path.isdir(path):
            print("Initializing directory {}".format(path))
            os.mkdir(path)

    return

