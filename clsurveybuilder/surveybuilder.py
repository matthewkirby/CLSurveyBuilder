"""Class to hold all of the grids and splines and performs the survey building"""
# Base python imports
import warnings
import numpy as np
from scipy.interpolate import interp2d

# This package
from . import observable
from . import masscatalog

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
    observables : list(Observable)
        List of observables of mock survey

    Notes
    -----
    - Implement A(z) rather than just A(z) = constant
    """

    def __init__(self, cfgfname, initialize_grids=True):
        """Create an instance of the SurveyBuilder class

        Parameters
        ----------
        cfgfname : str
            Name of the config file describing the survey realizations
        initialize_grids : bool (optional)
            If True, initializes 1D and 2D grids needed as identified by the config file
        """
        # Set the cosmology
        warnings.warn("Cosmology hard coded in SurveyBuilder", UserWarning)
        self.cosmo = {'flat': True, 'H0': 100., 'Om0': 0.300, 'Ob0': 0.046, 'sigma8': 0.80,
                      'ns': 0.972}
        setCosmology('survey cosmo', self.cosmo)

        # Manually set the redshift range of the survey
        warnings.warn("Redshift range for survey hard coded in SurveyBuilder", UserWarning)
        self.zmin = 0.2
        self.zmax = 0.65

        # Manually set survey area
        warnings.warn("Survey area constant and hard coded", UserWarning)
        self.solid_angle = 2500.*(np.pi**2/(180.**2))

        # Manually set the observables
        warnings.warn("SurveyBuilder is hard coded to SZObservable", UserWarning)
        self.observables = [observable.Observable(None)]

        if initialize_grids:
            self._initialize_grids()
            self.is_init = True
        else:
            self.is_init = False


    def _initialize_grids(self):
        r"""Initialize all of the required 1D and 2D grids for this survey realization.

        Notes
        -----
        - This function should identify which grids need to be initialized from the config file
        - Implements are m200m -> mdef spline
        - Inverse CDF for P(r_pearson) \propto 1-r^2 (TODO)
        - Halo mass function (if all reals are at fixed cosmology)? (TODO)
        - Inverse CDF for P(lamobs|lamtrue) from Costanzi19 (TODO)
        """
        self._init_mass_conversion_grids()


    def _init_mass_conversion_grids(self):
        """Build splines to convert M200m to the mass definitions that each observable is
        parameterized in. These splines are saved to the Observable object.
        TODO - Upper and lower mass ranges as optional params?
        """
        warnings.warn("Mass conversion spline is very coarse!")
        log10m200m = np.linspace(9, 17.5, 100, endpoint=False)
        # log10m200m = np.linspace(9, 17.5, 400, endpoint=False)
        m200m = np.power(10., log10m200m[np.where(log10m200m < 17.)])
        zgrid = np.linspace(self.zmin, self.zmax, 10)

        for obs in self.observables:
            if obs.mdef is '200m':
                continue

            # Build the 2d grid, spline it, set it to Observable()
            newmgrid = np.zeros((len(zgrid), len(m200m)))
            for i in range(len(zgrid)):
                c200m = concentration(m200m, '200m', zgrid[i])
                newmgrid[i,:] = changeMassDefinition(m200m, c200m, zgrid[i], '200m', obs.mdef)[0]
            obs.mass_conversion = interp2d(np.log(m200m), zgrid, np.log(newmgrid), kind='cubic',
                                           bounds_error=True)
        return


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












if __name__ == "__main__":
    cfgfname = ''
    x = SurveyBuilder(cfgfname, initialize_grids=True)
    x.build_survey('zetascatter', 1, save_m200=True, save_allmass=False, save_intr=False, save_obs=False)






