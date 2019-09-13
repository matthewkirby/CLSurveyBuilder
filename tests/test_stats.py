import sys, os
sys.path.append(os.environ['RICHESTRMCL'])
import configparser as cp
import numpy as np
from colossus.halo.mass_defs import changeMassDefinition
from colossus.cosmology.cosmology import setCosmology
from colossus.halo.concentration import concentration
from scipy.stats.stats import pearsonr
from scipy.interpolate import interp1d
from pandas import DataFrame

import mock_setup as ms
import draw_true_properties as dtp

import pytest



#@pytest.mark.skip(reason="Skipping for speed")
def test_draw_truth():
    # Load the fiducial model
    modelcfg = cp.ConfigParser()
    modelcfg.read('testing/mock_model0.ini')
    model = ms.build_model_dict(modelcfg)
    model['rich_pivot'] = np.log(modelcfg['Model'].getfloat('rich_pivot'))
    model['mgas_pivot'] = np.log(modelcfg['Model'].getfloat('mgas_pivot'))
    model['stdobs_lnmgas'] = modelcfg['Model'].getfloat('sig_mg_obs')

    # Set up the mass conversion
    cosparams = {'flat': True, 'H0': 100., 'Om0': 0.286, 'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.96}
    cosmo = setCosmology('myCosmo', cosparams)

    # Load the mass list
    n_objs = 10000000#0
    lnm200m_list = np.ones(n_objs)*np.log(2.3e14)
    lnm200m = lnm200m_list[0]
    lnm500c = np.log(changeMassDefinition(np.exp(lnm200m), 
                     concentration(np.exp(lnm200m), '200m', 0.2), 0.2, '200m', '500c')[0])
    lnm500c_list = np.ones(n_objs)*lnm500c
    catalog = DataFrame(np.array([lnm200m_list, lnm500c_list]).transpose(), columns=['lnm200m', 'lnm500c'])

    # Calculate truths
    lamtrue = 1.0 + model['lam0']*(np.exp(lnm200m-model['rich_pivot']))**model['alphalam']
    lamvar = (lamtrue*model['sigma0lam'])**2 + lamtrue
    mgtrue = model['mgas0']*(np.exp(lnm500c-model['mgas_pivot']))**model['alphamg']
    mgvar = (mgtrue*model['sigma0mg'])**2
    rtrue = model['r']

    # Draw the truths
    catalog = dtp.draw_true_richness(catalog, model)
    catalog = dtp.draw_true_mgas(catalog, model)

    # Print everything
    print("Mean Richness: {}, Truth: {}".format(np.mean(catalog['true_richness']), lamtrue))
    print("Richness STD: {}, Truth: {}".format(np.var(catalog['true_richness']), lamvar))
    print("Mean Mgas: {}, Truth: {}".format(np.mean(catalog['true_mgas']), mgtrue))
    print("Mgas STD: {}, Truth: {}".format(np.var(catalog['true_mgas']), mgvar))
    print("R: {}, Truth: {}".format(pearsonr(catalog['true_richness'], catalog['true_mgas'])[0], rtrue))

    # Make assertions
    np.testing.assert_almost_equal(np.mean(catalog['true_richness']), lamtrue, decimal=2)
    print("This test sometimes fails at the smaller n size")
    np.testing.assert_almost_equal(np.var(catalog['true_richness']), lamvar, decimal=1)
    np.testing.assert_almost_equal(np.mean(catalog['true_mgas']), mgtrue, decimal=2)
    np.testing.assert_almost_equal(np.var(catalog['true_mgas']), mgvar, decimal=1)
    np.testing.assert_almost_equal(pearsonr(catalog['true_richness'], catalog['true_mgas'])[0], rtrue, decimal=2)


@pytest.mark.skip(reason="This isn't actually testing anything..")
def test_obs_rich():
    n_objs = 10000000#0
    true_richness = np.ones(n_objs)*95.
    catalog = Table([true_richness.transpose()], names=(['true_richness']))
    invcdf_grid = ms.load_cdfs()
    catalog = dcp.draw_observed_richness(catalog, invcdf_grid)


@pytest.mark.skip(reason="no way of currently testing this")
def test_draw_obs():
    # Load the fiducial model
    model = cp.ConfigParser()
    model.read('testing/mock_model0.ini')
    model = dcp.build_model_dict(model)

    # Set up the mass conversion
    cosparams = {'flat': True, 'H0': 100., 'Om0': 0.286, 'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.96}
    cosmo = setCosmology('myCosmo', cosparams)

    # Load the mass list
    n_objs = 10000000#0
    lnm200m_list = np.ones(n_objs)*np.log(2.3e14)
    lnm200m = lnm200m_list[0]
    lnm500c = np.log(changeMassDefinition(np.exp(lnm200m), 
                     concentration(np.exp(lnm200m), '200m', 0.2), 0.2, '200m', '500c')[0])
    lnm500c_list = np.ones(n_objs)*lnm500c

    # Calculate truths
    lamtrue = model['lam0']*(np.exp(lnm200m-model['rich_pivot']))**model['alphalam']
    lamvar = (lamtrue*model['sigma0lam'])**2 + lamtrue
    mgtrue = model['mgas0']*(np.exp(lnm500c-model['mgas_pivot']))**model['alphamg']
    mgvar = (mgtrue*model['sigma0mg'])**2

    # Draw the truths
    true_richness, true_mgas, statistics = dcp.draw_truth(lnm200m_list, lnm500c_list, model)
    
    # Draw obs quantities
    varobs_lookup = dcp.load_varobs()
    observational_statistics = dcp.draw_observational_scatter(true_richness, model, 
                                                              statistics, varobs_lookup)
    obs_richness_adjustment = observational_statistics['obs_richness_adjustment']
    obs_mgas_adjustment = observational_statistics['obs_mgas_adjustment']
    obs_std_richness = observational_statistics['stdobs_richness']
    obs_std_mgas = observational_statistics['stdobs_mgas']

    obs_richness_scatter_truth = np.sqrt(interp1d(varobs_lookup[0], varobs_lookup[1], kind='cubic')(lamtrue))
    obs_mgas_scatter_truth = mgtrue*0.14

    np.testing.assert_almost_equal(np.mean(obs_richness_adjustment), 0.0, decimal=2)
    np.testing.assert_almost_equal(np.mean(obs_mgas_adjustment), 0.0, decimal=2)
    np.testing.assert_almost_equal(np.sqrt(np.var(obs_richness_adjustment)), obs_richness_scatter_truth, decimal=2)
    np.testing.assert_almost_equal(np.sqrt(np.var(obs_mgas_adjustment)), obs_mgas_scatter_truth, decimal=2)


if __name__ == "__main__":
    test_draw_truth()
