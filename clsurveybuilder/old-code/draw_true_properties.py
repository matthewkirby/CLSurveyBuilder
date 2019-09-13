"""
This module takes an mass catalog (200m) and draws true
richnesses and gas masses for each cluster.
"""
import numpy as np
from scipy.interpolate import UnivariateSpline, interp1d
import hmf_library as hmflib
import fit_properties as fp
import configparser as cp
from pandas import DataFrame
import sys
import matplotlib.pyplot as plt

import math
from time import time


def draw_true_mgas(catalog, model):
    """ Draw true richness for each cluster in the mass catalog

    Parameters
    ----------
    catalog : DataFrame
        Data frame with a row per cluster
    model : dict
        Dictionary of model parameters

    Returns
    -------
    Data frame with new columns for true richness
    """
    amp = model['mgas0']
    slope = model['alphamg']
    sigintr = model['sigma0mg']
    lnpivot = model['mgas_pivot']
    corr = model['r']

    mean = amp*np.exp((catalog['lnm500c']-lnpivot)*slope)
    var = np.power(mean*sigintr, 2)
    corr_mean = mean + np.sqrt(var/catalog['varint_richness'])*corr* \
                (catalog['true_richness']-catalog['mean_richness'])
    corr_var = (1.-corr**2)*var
    true_mgas = np.random.normal(corr_mean, np.sqrt(corr_var))

    catalog['mean_mgas'] = mean
    catalog['varint_mgas'] = var
    catalog['true_mgas'] = true_mgas

    return catalog


def draw_true_observable_from_file(catalog, model, filepath, rescale):
    corr = model['r']

    hinv = 1.0/0.7
    catalog['lnm500c_Msun'] = np.log(np.exp(catalog['lnm500c'])*hinv)

    catsize = float(len(catalog))
    catalog = catalog[catalog['lnm500c_Msun'] > np.log(1.e13)]
    catalog = catalog[catalog['lnm500c_Msun'] < np.log(1.e16)]
    cutsize = float(len(catalog))
    print("Cut: {} or {}\%".format(catsize-cutsize, (catsize/cutsize)-1.))

    # Set up the scaling relation lookup table
    mlookup, lxlookup = np.loadtxt(filepath, skiprows=1, usecols=[0,1], unpack=True)
    lnmlookup = np.log(mlookup)
    lnlxlookup = np.log(lxlookup)
    meanrel = interp1d(lnmlookup, lnlxlookup, kind='cubic')

    # Find the expected observed richness in mass bins
    cat, fitdata = fp.compute_expectation_bins(catalog, 'lnm200m', 'obs_richness', nbins=30)

    # Compute the means
    meanlnLx = meanrel(catalog['lnm500c_Msun'])
    varlnLx = model['sig0obs1']**2
    corr_meanlnLx = meanlnLx + np.sqrt(varlnLx/catalog['var_obsrichness'])*corr* \
                    (catalog['obs_richness']-catalog['expected_obs_richness'])
    corr_varlnLx = (1.-corr**2)*varlnLx
    true_lnLx = np.random.normal(corr_meanlnLx, np.sqrt(corr_varlnLx))
    
    # mean = np.exp(lnmean)/rescale
    # var = np.power(model['sig0obs1']*mean, 2)
    # corr_mean = mean + np.sqrt(var/catalog['varint_richness'])*corr* \
    #             (catalog['true_richness']-catalog['mean_richness'])
    # corr_var = (1.-corr**2)*var
    # true_Lx = np.random.normal(corr_mean, np.sqrt(corr_var))

    catalog['mean_truelnLx'] = meanlnLx
    catalog['var_truelnLx'] = varlnLx
    catalog['true_lnLx'] = true_lnLx

    return catalog



def draw_true_richness(catalog, model, parameterization='powerlaw', scattermodel='normal'):
    """ Draw true richness for each cluster in the mass catalog

    Parameters
    ----------
    catalog : DataFrame
        Data frame with a row per cluster
    model : dict
        Dictionary of model parameters
    parameterization : str, optional
        Use classic power law parameterization or Matteo's version
        using Mmin instead of lam0

    Returns
    -------
    Data frame with new columns for true richness
    """
    if parameterization == 'powerlaw':
        amp = model['lam0']
        slope = model['alphalam']
        sigintr = model['sigma0lam']
        lnpivot = model['rich_pivot']

        mean_richness = 1. + amp*np.exp((catalog['lnm200m']-lnpivot)*slope)
        scatter_intr = mean_richness*sigintr
        true_richness = np.random.poisson(mean_richness) + np.random.normal(0.0, scatter_intr)

        catalog['mean_richness'] = mean_richness
        catalog['varint_richness'] = mean_richness + scatter_intr**2
        catalog['true_richness'] = true_richness

    elif parameterization == 'matteo_for_tom':
        Mmin = model['mmin']
        Mpivot = model['mpivot']
        slope = model['alphalam']
        sigintr = model['sigma0lam']

        sat_richness = ((np.exp(catalog['lnm200m']) - Mmin)/(Mpivot))**slope
        scatter_intr = sat_richness*sigintr
        true_richness = 1.0 + np.random.poisson(sat_richness) + np.random.normal(0.0, scatter_intr)

        catalog['mean_richness'] = 1.0 + sat_richness
        catalog['varint_richness'] = sat_richness + scatter_intr**2
        catalog['true_richness'] = true_richness

    elif parameterization == 'matteo':
        Mmin = math.pow(10., model['logmmin'])
        M1 = math.pow(10., model['logm1'])
        slope = model['alphalam']
        sigintr = model['sig0lam']
        sigext = 0.33

        sat_richness = ((np.exp(catalog['lnm200m']) - Mmin)/(M1 - Mmin))**slope
        scatter_intr = sat_richness*sigintr

        # ===============================================
        # Gaussian draws
        # if scattermodel == 'normal':
        #     print("Normal lamtrue draws")
        #     mu = sat_richness + np.random.normal(0.0, scatter_intr)
        #     mu[mu<0.0] = 0.0
        #     catalog['true_richness'] = 0.9999 + np.random.poisson(mu)
        #     catalog['var_truerichness'] = sat_richness + scatter_intr**2
        #     catalog['var_obsrichness'] = sat_richness + (sigintr**2 + sigext**2)*sat_richness**2
        if scattermodel == 'normal':
            print("Normal lamtrue draws")
            catalog['true_richness'] = 1.0 + np.random.poisson(sat_richness) + np.random.normal(0.0, scatter_intr)
            catalog[catalog['true_richness'] < 0.0] = 0.0
            catalog['var_truerichness'] = sat_richness + scatter_intr**2
            catalog['var_obsrichness'] = sat_richness + (sigintr**2 + sigext**2)*sat_richness**2


        # ===============================================
        # Log normal
        # elif scattermodel == 'lognormal':
        #     print("Lognormal lamtrue draws")
        #     var = (1./sat_richness) + sigintr*sigintr
        #     expected_lnlamsat = np.log(sat_richness) - 0.5*var
        #     lnmu = expected_lnlamsat + np.random.normal(0.0, sigintr, size=len(expected_lnlamsat))
        #     mu = np.exp(lnmu)
        #     catalog['true_richness'] = 0.9999 + np.random.poisson(mu)
        #     catalog['var_truerichness'] = var
        #     catalog['var_obsrichness'] = var + sigext**2
        elif scattermodel == 'lognormal':
            print("Lognormal lamtrue draws")
            catalog['true_richness'] = 1.0 + np.random.poisson(sat_richness) + np.exp(np.random.normal(0.0, sigintr))
            catalog[catalog['true_richness'] < 0.0] = 0.0
            # catalog['var_truerichness'] = (1.0/sat_richness) + sigintr**2
            # catalog['var_obsrichness'] = (1.0/sat_richness) + sigintr**2 + sigext**2
            catalog['var_truerichness'] = sat_richness**2*((1.0/sat_richness) + sigintr**2)
            catalog['var_obsrichness'] = sat_richness**2*((1.0/sat_richness) + sigintr**2 + sigext**2)

        catalog['mean_truerichness'] = 1.0 + sat_richness


    else:
        raise TypeError("True richness parameterization {} not supported".format(parameterization))

    return catalog


# =================================================================================================
def draw_true_properties(cfgin, model, sinfo, file_names, splines, draw_mgas=True):
    """Draw the observed richness and mgas for a given mass catalog.

    Parameters
    ----------
    cfgin : ConfigParser(object)
        Inputs to the mcmc from the ini config file
    model : dict
        Dictionary of model parameters
    sinfo : SurveyInfo(object)
        Cosmology and other survey specific numbers
    file_names : dict
        Dictionary containing all of the output file names
    splines : dict
        Dictionary containing the splines for mass_conversion and varobs
    draw_mgas : Bool (True)
        Tells the function whether or not it should draw true mgas

    Outputs
    -------
    Dataframe with the true cluster properties
    """
    masscat = np.load(file_names['mass_catalog']+'.npy')
    catalog = DataFrame(masscat, columns=['lnm200m', 'z'])
    catalog = draw_true_richness(catalog, model)
    if draw_mgas:
        catalog['lnm500c'] = np.array([splines['mass_conversion'](row['lnm200m'], row['z']) for idx, row in catalog.iterrows()])
        catalog = draw_true_mgas(catalog, model)
    return catalog


# =================================================================================================
def find_richest_n(catalog, cfgin, file_names):
    """Given a list of halo observables, find the n halos that have the largest richhess

    Parameters
    ----------
    catalog : Dataframe
        Pandas dataframe containing the catalog info
    cfgin : ConfigParser(object)
        Inputs to the mcmc from the ini config file
    file_names : dict
        Dictionary containing all of the output file names

    Outputs
    -------
    The list of halo observable properties that got passed in for the richnest n halos
    """
    n_richest = cfgin['SurveyInfo'].getint('nclusters')
    catalog.sort_values('obs_richness', ascending=False, inplace=True)
    richest = catalog[:n_richest]
    richest = richest[['obs_richness', 'obs_richness_std', 'obs_mgas', 'obs_mgas_std', 'lnm200m']]

    # Make the output have reasonable length floats
    richest.to_csv(file_names['richest'], sep='\t', index_label='id', float_format='%.6e')

    return










