import sys, os
sys.path.append(os.environ['RICHESTRMCL'])
import numpy as np
import hmf_library as hmflib
import generate_mass_catalog as gmc
import configparser as cp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.integrate import quad
'''
@pytest.mark.skip(reason="takes too long and isnt the best test")
def test_mass_draws():
    return
    # Load the mcmc.ini file
    cfgin = cp.ConfigParser()
    cfgin.read('testing/mcmc.ini')
    sinfo = hmflib.SurveyInfo(cfgin)

    # Load a fiducial model
    model = cp.ConfigParser()
    model.read('testing/mock_model0.ini')
    deltaHMF = model['Model'].getfloat('dhmf')

    # Build a file names dictionary
    file_names = {"mass_catalog" : "testing/mass_cat.npy"}

    # Generate the mass catalog
    gmc.generate_mass_catalog(cfgin, model, file_names)
    lnmlist = np.load(file_names['mass_catalog'])

    # Load the halo mass function
    hmfpath = '../inputs/generateHMF/cosmosisOutput/mass_function/'
    lnmlist_grid, mhmf_grid = hmflib.read_hmf_files(hmfpath, cfgin, 1, 14.5)
    mhmf_interp = interp1d(lnmlist_grid, mhmf_grid*(1.+deltaHMF), kind='cubic')

    # Make mass bins
    dlnm = 0.01
    binedges = np.log(np.logspace(lnmlist[0]/np.log(10), lnmlist[-1]/np.log(10), num=int((lnmlist[-1]-lnmlist[0])/dlnm)))

    # Calculate the binned density
    binned_num, _ = np.histogram(lnmlist, binedges)

    # Calculate the expected binned density
    binned_num_truth = [sinfo.Vs*quad(mhmf_interp, binedges[i-1], binedges[i], epsabs=0, epsrel=1.e-8)[0]
                        for i in range(1, len(binedges))]

    # Calculate diff
    diff_in_bin = abs(binned_num-binned_num_truth)

    print diff_in_bin
    print 5*np.sqrt(binned_num_truth)

    # assert in each bin
    np.testing.assert_array_less(diff_in_bin, 5.*np.sqrt(binned_num_truth))
'''




if __name__ == "__main__":
    test_mass_draws()
