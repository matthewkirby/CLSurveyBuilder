"""Functions making general cosmology calculations"""
import numpy as np
from scipy.integrate import quad
import math

CLIGHT_KMS = 300000.

def hubble_parameter(z, cosmo):
    """Calculate the Hubble parameter at z for cosmology cosmo"""
    OmM, H0 = cosmo['Om0'], cosmo['H0']
    if cosmo['flat']:
        OmL = 1.0-OmM
        OmK = 0.0
    else:
        raise ValueError("Hubble parameter not defined for curved cosmology")
    return H0*math.sqrt(OmM*pow(1.+z, 3) + OmK*pow(1.+z, 2) + OmL)


def inv_hubble_parameter(z, cosmo):
    """Simple wrapper to calculate 1/H(z)"""
    return 1./hubble_parameter(z, cosmo)

def comoving_distance(z, cosmo):
    """Compute the comoving distance at z for cosmology cosmo"""
    return CLIGHT_KMS*quad(inv_hubble_parameter, 0.0, z, args=(cosmo), epsabs=1.49e-10,
                           epsrel=1.49e-10)[0]
