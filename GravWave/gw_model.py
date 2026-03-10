"""
gw_model.py

GW model data generator — to be called inside CoreModule in CoredataGenerator.py.

Takes sfr_popII (already computed by the reionization model) and the current
cosmological parameters, and computes the GW merger rate R(z_m) and the
luminosity distance grid DL(z_m).

These are stored in ctx so the GW likelihood in Likelihood.py can read them,
exactly mirroring how cl_TT, QHII etc. are stored by CoreModule.

Usage in CoredataGenerator.py
------------------------------
    from gw_model import GWModel

    gw_model = GWModel()   # instantiate once at module level

    # Inside CoreModule.__call__, after sfr_popII is computed:
    gw_model.compute(ctx, sfr_popII, H0, Om0)
    # This adds to ctx:
    #   'gw_R_merger'  — R(z_m) on z_grid_gw
    #   'gw_DL_grid'   — DL(z_m) on z_grid_gw
    #   'gw_z_grid'    — z_grid_gw itself
    #   'gw_dVc_grid'  — dVc/dz on z_grid_gw
"""

import numpy as np
import scipy
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM


# ---------------------------------------------------------------------------
# Reionization redshift grid — must match CoreModule Z array exactly
# ---------------------------------------------------------------------------
_zstart, _zend, _dz = 30.0, 0.0, 0.2
_Z_arraySize = int(round(abs((_zend - _zstart) / _dz)))
_Z_reion = np.array([
    ik * (-abs(_dz)) + _zstart for ik in range(_Z_arraySize + 1)
])
# descending: 30 → 0


class GWModel:
    """
    Computes GW model observables from the reionization model output.

    Parameters
    ----------
    z_min, z_max : float
        Redshift range for BBH merger rate integral.
    N_zgrid : int
        Number of points in GW redshift grid.
    t_delay : float
        BBH delay time in Gyr (delta-function, Eq 11).
    """

    def __init__(
        self,
        z_grid_gw,
        t_delay = None,
    ):
        self.t_delay = t_delay
        # GW merger redshift grid — fixed, cosmology-independent
        self.z_grid_gw = z_grid_gw


    def hubble_gyr(self, z, H0, ombh2, omch2):
         """H(z) in Gyr^-1 for flat LCDM."""
         h   = H0 / 100.0
         Om0 = (ombh2 + omch2) / h**2
         OmL = 1.0 - Om0
         E_z = np.sqrt(Om0 * (1 + z)**3 + OmL)
         return H0 * E_z * 1.022e-3   # Gyr^-1

    # ------------------------------------------------------------------
    # Main entry point — called inside CoreModule
    # ------------------------------------------------------------------

    def compute(self, H0, ombh2, omch2, Z_reion, sfr_popII, age_Gyr, D_L, dvc_dz):
        """
        Compute GW model outputs and store them in ctx.

        Parameters
        ----------
        ctx       : ChainContext — same context used by CoreModule
        sfr_popII : array, shape (N_Z_reion,)
            Pop II SFRD from reionization model on _Z_reion grid (30→0).
            This is Psi(z) driving the BBH merger rate via Eq 12.
        age_Gyr   :
        D_L       :
        dvc_dz    :
        """
        try:
            # --- Cosmological grids on z_grid_gw ---
            z_to_t_interp = interp1d(Z_reion, age_Gyr)           
            DL_interp   = interp1d(Z_reion, D_L, fill_value   = 'extrapolate')
            dVc_interp = interp1d(Z_reion, dvc_dz, fill_value   = 'extrapolate')
            t_to_z_interp = interp1d(age_Gyr, Z_reion, fill_value   = 'extrapolate')
            _sfr_log     = np.log10(np.maximum(sfr_popII, 1e-30))
            _sfrd_interp = interp1d(
                           Z_reion, _sfr_log,
                           kind         = 'linear',
                           bounds_error = False,
                           fill_value   = 'extrapolate')
            
            
            # Formation times for each merger redshift
            t_merger = z_to_t_interp(self.z_grid_gw)   # Gyr
            t_form   = t_merger - self.t_delay         

            valid_grid = t_form > 0
            zf_grid    = np.where(
                valid_grid,
                t_to_z_interp(np.maximum(t_form, 1e-3)),
                0.0
            )
            dzf_dtf = self.hubble_gyr(zf_grid, H0, ombh2, omch2) / (1.0 + zf_grid)/1e09

            # script_R(z_m) = Psi(z_f) * dz_f/dt_f  — Eq 11
            script_R_unnorm = np.where(valid_grid, 10**_sfrd_interp(zf_grid) * dzf_dtf, 0.0)
            # --- Normalise: script_R(z_m=4.822) = 20 Gpc^-3 yr^-1 ---
            script_R_at_4822 = interp1d(self.z_grid_gw, script_R_unnorm,
                                        bounds_error=False,
                                        fill_value='extrapolate')(4.822)
            norm             = 20.0 / script_R_at_4822
            
            # --- Observer-frame merger rate — Eq 10 ---
            DL_grid  = DL_interp(self.z_grid_gw)                          # D_L at merger z
            dVc_grid = dVc_interp(self.z_grid_gw) 
            R_merger = norm*(dVc_grid/1e09) / (1 + self.z_grid_gw) * script_R_unnorm #/yr/sr

            return R_merger, DL_grid, self.z_grid_gw
            
        except Exception as e:
            raise Exception(f"GWModel.compute failed: {e}")
            
