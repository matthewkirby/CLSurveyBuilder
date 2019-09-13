"""Class to fully describe an observable"""
import warnings

class Observable(object):
    """Fully describes an observable and can be referenced by SurveyBuilder

    Attributes
    ----------
    name : str
        Observable name used as header in outputs
    mscaling : str
        Mass scaling of the observable--mass relation. Chose an option from ['powerlaw']
        TODO - Implement 'hod', 'none'
    mdef : str
        Spherical overdensity mass definition for the mass scaling. Format is '200m' where the
        options are with respect to the mean 'm' or critical 'c' density. Unneeded/unused if
        mscaling is 'none'.
    zscaling : str
        Redshift scaling of the observable--mass relation Chose an option from ['Hz']
        TODO - Implement 'none', '1+z', 'Ez'
    intr_scatter_model : str
        Scatter model for intrinsic quantities given scaling relations. Chose an option from
        ['gauss', 'lognormal'].
        TODO - Implement 'none'
    obs_scatter_model : str
        Scatter model for observed quantities about intrinsic quantities. Chose an option from
        ['xizeta'].
        TODO - Implement 'C19', 'gauss', 'none'

    Notes
    -----
    - Should have a way to store constants for each scaling, pivot M and z for instance
    - Should have a way to store a prior on each model parameter
    - Should be loaded from a config file
    - Should know that I may append a mass conversion spline to it in SurveyBuilder()._init_m_conv
    """

    def __init__(self, cfg):
        """Create an instance of the Observable class

        Parameters
        ----------
        cfg :
            Config information for the survey realizations
        """
        warnings.warn("Observable class is hard coded for SZ observable", UserWarning)
        self.name = 'xi'
        self.mscaling = 'powerlaw'
        self.mdef = '500c'
        self.zscaling = 'Hz'
        self.intr_scatter_model = 'gauss'
        self.obs_scatter_model = 'xizeta'
