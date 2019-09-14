"""Class to hold all of the grids and splines and performs the survey building"""
# Base python imports
import warnings
from configparser import ConfigParser
import numpy as np
from scipy.interpolate import interp2d

# This package
try:
    from . import observable
    from . import masscatalog
except ImportError:
    import observable
    import masscatalog

# Colossus
from colossus.halo.mass_defs import changeMassDefinition
from colossus.cosmology.cosmology import setCosmology
from colossus.halo.concentration import concentration


class SurveyBuilder(object):
    """Class that initilizes the survey and holds precomputed quantities.

    The main object to generate surveys. It contains all precomputed quantities so that many
    realizations can be built without having to repeat slow computations.

    Attributes
    ----------
    outputpath : str
        Where do save the survey output
    surveyname : str
        Name of the survey. Output file will be {surveyname}{runid}.cat
    solid_angle : float
        Survey area in steradians
    zmin : float
        Minimum survey redshift
    zmax : float
        Maximum survey redshift
    cosmo : dict
        Cosmology parameters
    observables : list(Observable)
        List of observables of mock survey

    Notes
    -----
    - Implement A(z) rather than just A(z) = constant
    """

    def __init__(self, cfgfname):
        """Create an instance of the SurveyBuilder class

        Parameters
        ----------
        cfgfname : str
            Name of the config file describing the survey realizations
        """
        # Load the configparser object
        cfg = ConfigParser()
        cfg.read(cfgfname)

        # Set names and paths
        self.outputpath = cfg.get('General', 'surveypath')
        self.surveyname = cfg.get('General', 'surveyname')

        # Set survey properties
        self.solid_angle = cfg.getfloat('SurveyProps', 'area')*(np.pi**2/180.**2)
        self.zmin = cfg.getfloat('SurveyProps', 'zmin')
        self.zmax = cfg.getfloat('SurveyProps', 'zmax')

        # Load the cosmology model
        self._load_cosmology(cfg.get('General', 'cosmofile'))
        setCosmology('survey cosmo', self.cosmo)

        # Load the observables
        self.observables = []
        for i in range(1, cfg.getint('General', 'nobservables')+1):
            self.observables.append(observable.Observable(cfg, i, self, quick_cm_rel=True))


    def _load_cosmology(self, cosmofile):
        """Load the cosmology model and any scaling relation parameters

        Parameters
        ----------
        cosmofile : str
            Name of the cosmology file
        """
        cosmoin = ConfigParser()
        cosmoin.read(cosmofile)
        self.cosmo = {
                'flat': cosmoin.getboolean('Cosmology', 'flat'),
                'H0': cosmoin.getfloat('Cosmology', 'H0'),
                'Om0': cosmoin.getfloat('Cosmology', 'Omega_M'),
                'Ob0': cosmoin.getfloat('Cosmology', 'Omega_b'),
                'sigma8': cosmoin.getfloat('Cosmology', 'sigma8'),
                'ns': cosmoin.getfloat('Cosmology', 'ns')
        }
        return


    def _initialize_grids(self):
        r"""Initialize all of the required 1D and 2D grids for this survey realization.

        Notes
        -----
        - Inverse CDF for P(lamobs|lamtrue) from Costanzi19 (TODO)
        """
        pass


    def _init_realization(self, runi):
        """Setup that needs to be run once per realization, file names and model parameters

        Parameters
        ----------
        output : str
            Path to output directory for mocks
        runi : int
            Index indicating the ith run
        """
        return None, None
        

    def build_survey(self, catalogname, run_index, outpath='', save_m200=False, save_allmass=False,
                     save_intr=False, save_obs=False):
        """Build a mock cluster survey given with the desired outputs

        Parameters
        ----------
        catalogname : str
            File name of the catalog. Output will be in {catalogname}_{run_index}.extension
        run_index : int
            Index indicated the ith run. This will be appended to outputs for multiple
            realizations
        output : str (optional
            Additional prepended path to output directory for mocks
        save_m200 : bool
            Save cluster mass in m200 in the output table
        save_allmass : bool
            Save all masses used to compute scaling relations
        save_intr : bool
            Save intrinsic values of the scaling relation parameters in output
        save_obs : bool
            Save observed values of the observables in output
        """
        # fdict, model = self._init_realization(output, run_index)
        warnings.warn("Model parameters hard coded in build survey", UserWarning)
        model = {'s': 0.0, 'q': 1.0, 'Asz': 4.08, 'Bsz': 1.65, 'Csz': 0.64, 'sig0sz': 0.20}

        # Generate mass cat
        if save_m200 or save_allmass or save_intr or save_obs:
            print("Generating Mass Catalog")
            mock = masscatalog.generate_mass_catalog(self, model)
            pass

        # Generate intr cat
        if save_intr or save_obs:
            print("gen intr")
            pass

        # Generate obs cat
        if save_obs:
            print("gen obs")
            pass

        # Save the desired outputs

        return mock






