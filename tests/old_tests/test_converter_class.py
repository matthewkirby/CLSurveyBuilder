import sys, os
sys.path.append(os.environ['RICHESTRMCL'])
import numpy as np
import draw_observed_properties as dop
import pytest

@pytest.mark.skip(reason="Skipping for speed")
def test_converter():
    x = dop.Converter()

    # Testing that the object was made
    assert x.is_setup == False
    assert x.icdf_is_setup == False

    # Test setting up model
    x.setup_splines()
    assert x.is_setup == True
    assert x.icdf_is_setup == False

    # Test the probability distribution to 5%
    ltrue_test = 95
    lobs_grid, probinput = np.loadtxt('example_proj_model/lamtrue{}_plamobs.txt'.format(ltrue_test), unpack=True)
    probclass = x.p_lambda_obs_true(ltrue_test, lobs_grid, 0.20)
    np.testing.assert_allclose(probclass, probinput, rtol=5.e-2)

    # Test the CDF
    x.setup_iCDF_grid()
    assert False
