"""Class to fully describe an observable"""
import numpy as np
from scipy.interpolate import interp2d
from colossus.halo.mass_defs import changeMassDefinition
from colossus.cosmology.cosmology import setCosmology
from colossus.halo.concentration import concentration


class Observable(object):
    """Fully describes an observable and can be referenced by SurveyBuilder

    Attributes
    ----------
    name : str
        Observable name used as header in outputs
    scaling : str
        Model of the observable--mass relation.
        Options: powerlaw_mz, lambda_hod
    mdef : str
        Spherical overdensity mass definition for the mass scaling.
    intr_scatter_model : str
        Scatter model for intrinsic quantities given scaling relations.
        Options: gaussian, lognormal, none
    obs_scatter_model : str
        Scatter model for observed quantities about intrinsic quantities.
        Options: sz, costanzi19
    """

    def __init__(self, cfg, obsind, builder, quick_cm_rel=False):
        """Create an instance of the Observable class

        Parameters
        ----------
        cfg : ConfigParser
            Config information for the survey realizations
        obsind : int
            Index of the observable. Loads from cfg section Observable{obsind}
        builder : SurveyBuilder
            The SurveyBuilder object calling this constructor. Used for survey properties
            for the mass conversion spline
        quick_cm_rel : bool
            Optional flag to make a coarse grid for the c-M relation that builds much faster
        """
        secname = 'Observable{}'.format(int(obsind))
        self.name = cfg.get(secname, 'name')
        self.scaling = cfg.get(secname, 'scaling')
        self.mdef = cfg.get(secname, 'mdef')
        self.intr_scatter_model = cfg.get(secname, 'intr_scatter_model')
        self.obs_scatter_model = cfg.get(secname, 'obs_scatter_model')

        self._config_validation(obsind)

        if self.mdef is '200m':
            self.mass_conversion = None
        else:
            self.mass_conversion = self._spline_mass_conversion(builder, quick_cm_rel)


    def _spline_mass_conversion(self, builder, quick):
        """2D Spline of the conversion from M200m to mdef"""
        if quick:
            npoints = 100
        else:
            npoints = 400
        log10m200m = np.linspace(9, 17.5, npoints, endpoint=False)
        m200m = np.power(10., log10m200m[np.where(log10m200m < 17.)])
        zgrid = np.linspace(builder.zmin, builder.zmax, 10)

        newmgrid = np.zeros((len(zgrid), len(m200m)))
        for i in range(len(zgrid)):
            c200m = concentration(m200m, '200m', zgrid[i])
            newmgrid[i,:] = changeMassDefinition(m200m, c200m, zgrid[i], '200m', self.mdef)[0]
        mass_conversion = interp2d(np.log(m200m), zgrid, np.log(newmgrid), kind='cubic',
                                   bounds_error=True)
        return mass_conversion


    def __repr__(self):
        outstr1 = '{}\n\tScaling: {}\n\tmdef: {}\n'.format(self.name, self.scaling, self.mdef)
        outstr2 = '\tintrinsic scatter model: {}\n'.format(self.intr_scatter_model)
        outstr3 = '\tobservational scatter model: {}\n'.format(self.obs_scatter_model)
        return outstr1+outstr2+outstr3


    def _config_validation(self, obsind):
        scaling = ['powerlaw_mz', 'lambda_hod']
        if self.scaling not in scaling:
            raise ValueError("{} scaling of observable {} invalid".format(self.scaling, obsind))
        if self.mdef[-1] not in ['m', 'c']:
            raise ValueError("{} not a valid SO mass definition".format(self.mdef))
        intrmodels = ['gaussian', 'lognormal', 'none']
        if self.intr_scatter_model not in intrmodels:
            raise ValueError("Invalid intrinsic scatter model for observable{}".format(obsind))
        obsmodels = ['sz', 'costanzi19']
        if self.obs_scatter_model not in obsmodels:
            raise ValueError("Invalid obs scatter model for observable{}".format(obsind))
        return
