from scipy.integrate import quad
import matplotlib.pyplot as plt
import hmf_library as hmflib
import numpy as np
import sys, os
from pandas import DataFrame
import pandas as pd
from scipy.stats import poisson

from progress.bar import Bar
from time import time

from colossus.cosmology.cosmology import setCosmology                                               
from colossus.lss.mass_function import massFunction   


def draw_clusters(grid):
    for idx, row in grid.iterrows():
        lnmlist = np.random.uniform(row['lnmlo'], row['lnmhi'], int(row['N']))
        zlist = np.random.uniform(row['zlo'], row['zhi'], int(row['N']))

        try:
            catalog = np.vstack((catalog, np.array([lnmlist, zlist]).transpose()))
        except NameError:
            catalog = np.array([lnmlist, zlist]).transpose()

    return DataFrame(catalog, columns=['lnm200m', 'z'])


def draw_nbin(grid):
    grid[grid['<N>'] < 0.0] = 0
    grid['N'] = np.random.poisson(grid['<N>'].values)
    return


def poisson_lnm_integrand(lnm, z, sinfo, model):
    shmf, qhmf = model['shmf'], model['qhmf']
    log10Mstar = 13.8 # hinv Mphc
    return massFunction(np.exp(lnm), z, q_in='M', q_out='dndlnM', mdef='200m', model='tinker08')*(shmf*(lnm/np.log(10.) - log10Mstar) + qhmf)


def poisson_z_integrand(z, row, sinfo, model):
    dV = hmflib.volumeIntegrand(z, sinfo)
    lnms = np.array([row['lnmlo'], (row['lnmlo']+row['lnmhi'])*0.5, row['lnmhi'] ])
    Ilnms = poisson_lnm_integrand(lnms, z, sinfo, model)
    #mass_integral = quad(poisson_lnm_integrand, row['lnmlo'], row['lnmhi'], args=(z, sinfo),
    #                     epsabs=1e-4)[0]
    #return dV*mass_integral
    return dV*np.trapz(Ilnms, lnms)


def compute_poisson_mean(grid, sinfo, model):
    t0 = time()
    bar = Bar('Computing <n>', max=len(grid))
    for idx, row in grid.iterrows():
        bar.next()
        row['<n>'] = quad(poisson_z_integrand, row['zlo'], row['zhi'], args=(row, sinfo, model),
                          epsabs=1e-4)[0]
    t1 = time()
    bar.finish()
    print("Computing the mean took {} seconds".format(t1-t0))
    return


def build_grids(lnmlist, cfgin):
    # print("USING COARSER GRID IN MASS CATALOG")
    dlog10m, dz = 0.04, 0.02
    # dlog10m, dz = 0.008, 0.002
    # dlog10m, dz = 0.004, 0.002
    log10mlist = lnmlist/np.log(10)
    mgrid = np.logspace(log10mlist[0], log10mlist[-1], num=int((log10mlist[-1]-log10mlist[0])/dlog10m))
    lnmgrid = np.log(mgrid)
    zgrid = np.arange(cfgin['SurveyInfo'].getfloat('zlo'), cfgin['SurveyInfo'].getfloat('zhi')+.000001, dz)

    # Make combinations and save in dataframe
    ndgrid = []
    for i in range(1, len(lnmgrid)):
        for j in range(1, len(zgrid)):
            ndgrid.append(np.array([lnmgrid[i-1], lnmgrid[i], zgrid[j-1], zgrid[j]]))
    ndgrid = np.array(ndgrid)

    # Add the empty columns and convert to data frame
    ndgrid = np.hstack((ndgrid, -1.*np.ones((ndgrid.shape[0], 3))))
    fullgrid = DataFrame(ndgrid, columns=['lnmlo', 'lnmhi', 'zlo', 'zhi', '<n>', '<N>', 'N'])

    return fullgrid


"""
REFACTOR COMMENT:
-----------------
>Once I make fullgrid and populate it with <N> and N, just remove any rows with N=0
rather than looping over them. This should give a pretty big speed up.
>If the grid already exists, the mass cut should cut the grid. If the mass cut is
lower than the grid minimum, give an error!
"""
def generate_mass_catalog(cfgin, model, sinfo, file_names, log10mcut=14., save=True):
    '''Generate a mass catalog given a survey volumne'''
    lnmlist, _ = hmflib.load_tinker_hmf(0.2, sinfo, -1., 17.)
    cosmology = {'flat': True, 'H0': sinfo.H0, 'Om0': sinfo.OmM, 'Ob0': sinfo.Omb,
                 'sigma8': sinfo.sigma8, 'ns': sinfo.ns}
    setCosmology('hmf_cosmo', cosmology)


    # =========================================================================================
    # Grid in M and z and permute in dataframe. Save it to a file
    # try:
    #     pklpath = os.environ['HMFGRIDPATH']
    # except KeyError:
    #     pklpath = 'grids'

    # zmin, zmax = cfgin['SurveyInfo']['zlo'], cfgin['SurveyInfo']['zhi']
    # pklname = 'precomputed_poisson_mean_{}_zlo{}_zhi{}_100s{:.3f}_q{:.3f}.pkl'.format(cfgin['SurveyInfo']['cosmology'], zmin, zmax, 100.*model['shmf'], model['qhmf'])
    # try:
    #     fullgrid = pd.read_pickle(os.path.join(pklpath,pklname))
    #     print("Loaded hmf grid")
    # except IOError:
    #     print("Grid not found...Computing grid with log10cut = {}".format(log10mcut))
    #     fullgrid = build_grids(np.array(lnmlist)[np.where(lnmlist > log10mcut*np.log(10))], cfgin)
    #     compute_poisson_mean(fullgrid, sinfo)
    #     fullgrid.to_pickle(os.path.join(pklpath,pklname))
    fullgrid = build_grids(np.array(lnmlist)[np.where(lnmlist > log10mcut*np.log(10))], cfgin)
    compute_poisson_mean(fullgrid, sinfo, model)


    # =========================================================================================
    # Apply delta_hmf and draw the number of objects per bin
    fullgrid['<N>'] = fullgrid['<n>']*sinfo.solid_angle
    draw_nbin(fullgrid)


    # =========================================================================================
    # Now build the catalog
    catalog = draw_clusters(fullgrid)
    print("Mean cluster redshift = {}".format(np.mean(catalog['z'])))

    # =================================================
    # Write out mass catalog
    if save:
        np.save(file_names['mass_catalog'], catalog.values)
        return
    else:
        return catalog
