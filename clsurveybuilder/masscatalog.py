"""Build a mass catalog from the halo mass function"""
import warnings
from time import time
from progress.bar import Bar
from scipy.integrate import quad
import numpy as np
from pandas import DataFrame


from colossus.lss.mass_function import massFunction

# Imports from this package
try:
    from . import cosmology_funcs as cfuncs
except ImportError:
    import cosmology_funcs as cfuncs


def _draw_clusters(grid):
    """Draw clusters in mass and redshift bins"""
    started = False
    for _, row in grid.iterrows():
        lnmlist = np.random.uniform(row['lnmlo'], row['lnmhi'], int(row['N']))
        zlist = np.random.uniform(row['zlo'], row['zhi'], int(row['N']))

        if not started:
            catalog = np.array([lnmlist, zlist]).transpose()
            started = True
        else:
            catalog = np.vstack((catalog, np.array([lnmlist, zlist]).transpose()))

    return DataFrame(catalog, columns=['lnm200m', 'z'])


def _n_given_z(lnm, redshift, model, cosmo):
    """Integrand over mass at fixed z for <n>"""
    if abs(cosmo['H0'] - 100.) > 0.1:
        raise ValueError("Mass Function to build mocks is in units of hinv Msun")

    # massfunc_options = {'q_in':'M', 'q_out':'dndlnM', 'mdef':'200m', 'model':'tinker08'}
    #
    shmf, qhmf = model['s'], model['q']
    nuisance = (shmf*(lnm/np.log(10.) - 13.8) + qhmf) # Pivot mass in Msun hinv
    return massFunction(np.exp(lnm), redshift, q_in='M', q_out='dndlnM',
                        mdef='200m', model='tinker08')*nuisance
    # return massFunction(np.exp(lnm), redshift, *massfunc_options)*nuisance


def _n_given_massz(z, row, model, builder):
    """Integrand over redshift of <n>"""
    chi_comov = cfuncs.comoving_distance(z, builder.cosmo)
    dvdzda = cfuncs.CLIGHT_KMS*cfuncs.inv_hubble_parameter(z, builder.cosmo)*chi_comov*chi_comov
    lnms = np.array([row['lnmlo'], (row['lnmlo']+row['lnmhi'])*0.5, row['lnmhi']])
    integrated_vals_lnms = _n_given_z(lnms, z, model, builder.cosmo)
    return dvdzda*np.trapz(integrated_vals_lnms, lnms)


def _compute_expected_density(grid, model, builder):
    """Compute expected halo density in each mass and redshift bin

    This is simple the integral of the mass function over the bin

    Parameters
    ----------
    grid : DataFrame
        A table of lnM and z for which to draw masses
    model : dict
        Model parameters to generate mock
    builder : SurveyBuilder
        Contains all of the survey information and splines
    """
    timestart = time()
    progbar = Bar('Computing <n>', max=len(grid))
    for _, row in grid.iterrows():
        progbar.next()
        row['<n>'] = quad(_n_given_massz, row['zlo'], row['zhi'], args=(row, model, builder),
                          epsabs=1e-4)[0]
    timeend = time()
    progbar.finish()
    print("Computing the mean took {} seconds".format(timeend-timestart))
    return


def _build_grids(lnmlist, builder):
    """Using the builder object and provided ln(M)'s, make a grid of _________

    Parameters
    ----------
    lnmlist : np.ndarray
        Natural log of the halo mass in 200m
    builder : SurveyBuilder
        Contains all of the survey information and splines
    """
    # Set the grid deltas for both log10mass and redshift
    warnings.warn("Using coarse grid for hmf calculation", UserWarning)
    deltalog10m = 0.1 # Historically 0.04, 0.008, 0.004
    deltaz = (builder.zmax-builder.zmin)/20. # Historically 0.02 or 0.002

    # Build the grid
    log10mmin, log10mmax = np.min(lnmlist)/np.log(10.), np.max(lnmlist)/np.log(10.)
    lnmgrid = np.log(np.logspace(log10mmin, log10mmax, num=int((log10mmax-log10mmin)/deltalog10m)))
    zgrid = np.arange(builder.zmin, builder.zmax+.000001, deltaz)

    # Make combinations - I think there is a numpy function to do this?
    ndgrid = []
    for i in range(1, len(lnmgrid)):
        for j in range(1, len(zgrid)):
            ndgrid.append(np.array([lnmgrid[i-1], lnmgrid[i], zgrid[j-1], zgrid[j]]))
    ndgrid = np.array(ndgrid)

    # Add the empty cols and save to data frame
    ndgrid = np.hstack((ndgrid, -1.*np.ones((ndgrid.shape[0], 3))))
    fullgrid = DataFrame(ndgrid, columns=['lnmlo', 'lnmhi', 'zlo', 'zhi', '<n>', '<N>', 'N'])

    return fullgrid


def generate_mass_catalog(builder, model, log10mcut=14.):
    """Generate a mass catalog given a survey volume

    Parameters
    ----------
    builder : SurveyBuilder
        Contains all of the survey information and splines
    model : dict
        Model parameters to generate mock
    log10mcut : float
        Lower mass cut for survey in 200m. This value should be chose so as to not
        impact the richest X clusters.

    Notes
    -----
    >Once I make fullgrid and populate it with <N> and N, just remove any rows with N=0
    rather than looping over them. This should give a pretty big speed up.
    >If the grid already exists, the mass cut should cut the grid. If the mass cut is
    lower than the grid minimum, give an error!
    """

    # Load the mass list and make cuts
    warnings.warn("Mass grid in gen mass catalog is pretty coarse", UserWarning)
    log10m200m = np.linspace(9, 17.5, 100, endpoint=False) # 400 bins?
    log10m200m = log10m200m[np.where((log10m200m > log10mcut) & (log10m200m < 17.))]
    lnmlist = np.log(np.power(10., log10m200m))

    # MK 9/12/2019 I don't know what this code is
    # =========================================================================================
    # Grid in M and z and permute in dataframe. Save it to a file
    # try:
    #     pklpath = os.environ['HMFGRIDPATH']
    # except KeyError:
    #     pklpath = 'grids'

    # zmin, zmax = cfgin['SurveyInfo']['zlo'], cfgin['SurveyInfo']['zhi']
    # pklname = 'precomputed_poisson_mean_{}_zlo{}_zhi{}_100s{:.3f}_q{:.3f}.pkl'.format(
    #             cfgin['SurveyInfo']['cosmology'], zmin, zmax, 100.*model['shmf'], model['qhmf'])
    # try:
    #     fullgrid = pd.read_pickle(os.path.join(pklpath,pklname))
    #     print("Loaded hmf grid")
    # except IOError:
    #     print("Grid not found...Computing grid with log10cut = {}".format(log10mcut))
    #     fullgrid = _build_grids(np.array(lnmlist)[np.where(lnmlist > log10mcut*np.log(10))],
    #                             cfgin)
    #     compute_poisson_mean(fullgrid, sinfo)
    #     fullgrid.to_pickle(os.path.join(pklpath,pklname))
    # =========================================================================================
    # MK 9/12/2019 End confusion

    # Compute the the hmf at each grid point
    fullgrid = _build_grids(lnmlist, builder)
    _compute_expected_density(fullgrid, model, builder)

    # =========================================================================================
    # Compute <N> (fold this into the previous step)
    fullgrid['<N>'] = fullgrid['<n>']*builder.solid_angle

    # =========================================================================================
    # Draw N per bin with some optimizations
    fullgrid['N'] = np.random.poisson(fullgrid['<N>'].values)

    # =========================================================================================
    # Now build the catalog
    catalog = _draw_clusters(fullgrid)
    print("Mean cluster redshift = {}".format(np.mean(catalog['z'])))

    # =================================================
    # Write out mass catalog
    return catalog
