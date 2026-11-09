from cobaya.likelihood import Likelihood as CobayaLikelihood
from typing import Optional
import numpy as np

from Likelihood import Likelihood as OldLikelihood
from CoredataGenerator import CoreModule


class ReionLikelihood(CobayaLikelihood):

    
    # These values are read from input.yaml
    include_planck: bool = True
    include_gw: bool = True
    include_21cm: bool = True

    gamma_datafile: Optional[str] = None
    lyman_limit_datafile: Optional[str] = None
    T_b_Datafile: Optional[str] = None

    input_params = [
        "hubble", "ombh2", "omch2", "As", "ns",
        "f_zero", "alpha_lo", "alpha_hi", "alpha_z",
        "log_esc_II", "lambda0", "log_f_X","log_f_alpha"
    ]

    output_params = ["Q_HII_at_z5p8", "tau"]
    
    def initialize(self):

        # Instantiate your original likelihood
        self.like = OldLikelihood(
            CoreModule,
            include_planck=self.include_planck,
            include_gw=self.include_gw,
            include_21cm=self.include_21cm,
        )

        # Initialise observational data
        self.like.setup(
            self.gamma_datafile,
            self.lyman_limit_datafile,
            self.T_b_Datafile,
        )

        self._derived = {}

    def logp(
        self,
        hubble,
        ombh2,
        omch2,
        As,
        ns,
        f_zero,
        alpha_lo,
        alpha_hi,
        alpha_z,
        log_esc_II,
        lambda0,
        log_f_X,
        log_f_alpha,
        _derived= None
    ):

        # Parameter vector in exactly the same order as emcee
        p = np.array([
            hubble,
            ombh2,
            omch2,
            As,
            ns,
            f_zero,
            alpha_lo,
            alpha_hi,
            alpha_z,
            log_esc_II,
            lambda0,
            log_f_X,
            log_f_alpha
        ])

        # Build ChainContext exactly as before
        ctx = self.like.createChainContext(p)

        # Run the expensive model
        ok = self.like.Core(ctx)

        if not ok:
            if _derived is not None:
                _derived["Q_HII_at_z5p8"] = np.nan
                _derived["tau"] = np.nan
            return -np.inf

        # Compute likelihood
        loglike, tau, Q_HII_at_z5p8 = self.like.ComputeLikelihood(ctx)
        
        if _derived is not None:
            _derived["Q_HII_at_z5p8"] = Q_HII_at_z5p8
            _derived["tau"] = tau
        
        return loglike



