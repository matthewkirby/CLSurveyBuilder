import numpy as np
import sys
import configparser as cp
import hmf_library as hmflib
from astropy.table import Table
from scipy.interpolate import interp2d, interp1d, griddata
import time
import matplotlib.pyplot as plt
import mock_setup as ms
from scipy.integrate import quad

# =================================================================================================
def compute_RM_given_MR(cfgin, sig0lam, amp, slope):
    cosmology = cfgin['SurveyInfo']['cosmology']
    table = np.loadtxt('../priors/lookup_table_v2/tables/final_lookup_{}.dat'.format(cosmology))
    rows = Table(table, names=('Alam', 'Blam', 'Slam', 'amp', 'slope'))
    subtable = rows[np.where(abs(rows['Slam'] - sig0lam) < 0.001)]

    points = [[r['amp'], r['slope']] for r in subtable]
    Alam = griddata(points, subtable['Alam'], (amp, slope), method='cubic')
    Blam = griddata(points, subtable['Blam'], (amp, slope), method='cubic')

    return Alam, Blam


# =================================================================================================
def draw_simet(priors):
    amp = np.random.normal(priors.amp_mr_rel, priors.amp_uncert_mr_rel)
    slope = np.random.normal(priors.slope_mr_rel, priors.slope_uncert_mr_rel)
    return amp, slope


# =================================================================================================
def rm_relation_from_prior(priors, sig0lam, cfgin):
    opts = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    sig0lamgrid = opts[np.argmin(abs(opts-sig0lam))]

    count = 0
    while True:
        count += 1
        amp, slope = draw_simet(priors)
        Alam, Blam = compute_RM_given_MR(cfgin, sig0lamgrid, amp, slope)
        if 30. < Alam < 90. and 0.5 < Blam < 1.:
            return Alam, Blam
        if count == 50:
            print("No parameters found in 50 tries.")
            sys.exit(1)
        print("Out of Bounds, redrawing")


# =================================================================================================
def rm_relation_from_flat():
    Alam = np.random.uniform(30., 90.)
    Blam = np.random.uniform(0.5, 1.0)
    return Alam, Blam


# =================================================================================================
def draw_mf_nuisance_params():
    mus = [0.037, 1.008]
    C = [[0.00019, 0.00024], [0.00024, 0.00038]]
    return np.random.multivariate_normal(mus, C)


# =================================================================================================
def draw_corr_coefficient(splines):
    """Assuming P(r) \propto 1 - r^2"""
    return splines['P(r)'](np.random.random())


# =================================================================================================
def draw_model(cfgin, sinfo, file_names, splines):
    '''Draw a fiducial model from priors

    Parameters
    ----------
    cfgin : ConfigParser(object)
        Inputs to the mcmc from the ini config file
    sinfo : SurveyInfo(object)
        Cosmology and other survey parameters
    file_names : dict
        Dictionary containing all of the output file names
    splines : dict
        Dictionary of spline lookup tables - Used for P(r)

    Returns
    -------
    A dictionary of model parameters
    '''
    model = cp.ConfigParser()
    model.add_section('Model')
    priors = hmflib.PriorContainer(cfgin, sinfo)
    priortype = cfgin['MCMC']['prior']

    # Draw the rest of the parameters
    model.set('Model', 'r', str(draw_corr_coefficient(splines)))

    # Set the Mgas-Mass relation
    model.set('Model', 'mgas0', str(np.random.normal(priors.mg0, priors.mg0_uncert)))
    model.set('Model', 'alphamg', str(np.random.normal(priors.alphamg, priors.alphamg_uncert)))
    model.set('Model', 's0mg', str(np.random.uniform(0.03, 0.11)))

    # Set the Richness-Mass relation
    fixed_vintr = cfgin['General'].getfloat('fixed_vintr')
    if fixed_vintr > 0.0:
        model.set('Model', 's0lam', str(np.sqrt(fixed_vintr)))
    else:
        model.set('Model', 's0lam', str(np.sqrt(np.random.uniform(0.01, 0.25))))
    # Alam, Blam = rm_relation_from_prior(priors, model['Model'].getfloat('s0lam'), cfgin)
    Alam, Blam = rm_relation_from_flat()
    model.set('Model', 'lam0', str(Alam))
    model.set('Model', 'alphalam', str(Blam))

    s, q = draw_mf_nuisance_params()
    model.set('Model', 'shmf', str(s))
    model.set('Model', 'qhmf', str(q))

    # Write the model to a file
    with open(file_names['model'], 'w') as modelout:
        model.write(modelout)

    mdict = ms.build_model_dict(model)

    return mdict


def print_model(model):
    print("r = {} \nDhmf = {} \nAlam = {} \nBlam = {} \nSlam = {} \nAmg = {} \nBmg = {} \nSMg = {}"
          .format(model['Model']['r'], model['Model']['Dhmf'], model['Model']['lam0'], 
                  model['Model']['alphalam'], model['Model']['s0lam'], model['Model']['mgas0'],
                  model['Model']['alphamg'], model['Model']['s0mg']))
    return


if __name__ == "__main__":
    cfgin = cp.ConfigParser()
    cfgin.read('../mcmc.ini')

    table = load_lookup_table(0.1, cfgin)
    plt.plot(table['amp'], table['slope'], 'o', label='0.1')
    table = load_lookup_table(0.225, cfgin)
    plt.plot(table['amp'], table['slope'], 'o', label='0.225')
    table = load_lookup_table(0.35, cfgin)
    plt.plot(table['amp'], table['slope'], 'o', label='0.35')
    table = load_lookup_table(0.475, cfgin)
    plt.plot(table['amp'], table['slope'], 'o', label='0.475')
    table = load_lookup_table(0.6, cfgin)
    plt.plot(table['amp'], table['slope'], 'o', label='0.6')
    plt.legend()
    plt.show()
    sys.exit(1)

    file_names = {'model' : 'blahblah.dat'}

    model = draw_model(cfgin, file_names)




