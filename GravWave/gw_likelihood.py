"""
gw_likelihood.py

GW hyperlikelihood — to be called inside Likelihood.ComputeLikelihood().

Reads GW model outputs stored in ctx by GWModel (via CoredataGenerator):
    ctx.get('gw_R_merger') — merger rate R(z_m) on z_grid_gw
    ctx.get('gw_DL_grid')  — luminosity distance DL(z_m) on z_grid_gw
    ctx.get('gw_z_grid')   — redshift grid

Computes Eq 18:
    ln L_GW = -zeta + sum_i ln L_i

Usage in Likelihood.py
-----------------------
    from gw_likelihood import GWLikelihood

    # In setup():
    self.gw_like = GWLikelihood()
    self.gw_like.setup()

    # In ComputeLikelihood():
    loglike_gw   = self.gw_like.ComputeLikelihood(ctx)
    loglike_tot += loglike_gw
"""

import numpy as np
import pickle
from scipy.integrate import simpson
import matplotlib.pyplot as plt



class GWLikelihood:
    """
    GW hyperlikelihood — Eq 18.

    Parameters
    ----------
    catalogue_path : str
    T_obs          : float — observation time in years
    """

    def __init__(
        self,
        #change path
        catalogue_path = "/Users/users/achatterjee/Data/GW_with_CosmoReionMC/GravWave/CE_mock_catalog.pkl",
        T_obs          = 70.0 / 365.25,
    ):
        self.catalogue_path = catalogue_path
        self.T_obs          = T_obs
        self._is_setup      = False

    # ------------------------------------------------------------------
    # Setup — load catalogue once
    # ------------------------------------------------------------------

    def setup(self, fiducial_H0=67.7, fiducial_Om0=0.31):
        """
        Load GW catalogue and precompute fiducial sigma_DL rescaling.

        Parameters
        ----------
        fiducial_H0, fiducial_Om0 : float
            Cosmology used when generating the catalogue.
        """
        with open(self.catalogue_path, "rb") as f:
            catalogue = pickle.load(f)

        self.N_events     = len(catalogue)
        self.DL_obs       = np.array([
            ev["injection"]["luminosity_distance"] for ev in catalogue
        ])                                                     # Mpc, shape (N,)
        self.sigma_DL_fid = np.array([ev["sigma_DL"] for ev in catalogue])
        self.z_true       = np.array([ev["z_true"] for ev in catalogue])
        self.R_obs        = np.array([ev["R_true"] for ev in catalogue])

        self._is_setup = True

    # ------------------------------------------------------------------
    # Likelihood — reads model outputs from ctx
    # ------------------------------------------------------------------

    def ComputeLikelihood(self, ctx):
        """
        Compute ln L_GW from ctx populated by GWModel in CoredataGenerator.

        Parameters
        ----------
        ctx : ChainContext

        Returns
        -------
        float : ln L_GW, or -inf if invalid
        """
        if not self._is_setup:
            raise RuntimeError("Call GWLikelihood.setup() before ComputeLikelihood().")

        # --- Read model outputs obtained by GWModel in each MCMC step---
        R_grid  = ctx.get('gw_R_merger')
        DL_grid = ctx.get('gw_DL_grid')
        z_grid  = ctx.get('gw_z_grid')
        
        if R_grid is None or DL_grid is None:
            return -np.inf

        if np.any(~np.isfinite(R_grid)) or np.any(R_grid < 0):
            return -np.inf
        norm = simpson(R_grid, z_grid)
        zeta = self.T_obs * norm #should we multiply with 4*pi
        #R_tilde = R_grid/norm
        
        if zeta <= 0 or not np.isfinite(zeta):
            return -np.inf
        
        # --- Rescale sigma_DL for current cosmology ---
        DL_model_at_ztrue = np.interp(self.z_true, z_grid, DL_grid)
        sigma_DL          = self.sigma_DL_fid * (
            DL_model_at_ztrue / self.DL_obs
        ) 
        
        #added this part to calculate chi^2 for merger rate
        R_model_at_ztrue        =  np.interp(self.z_true, z_grid, R_grid)
        sigma_R                 = np.sqrt(self.R_obs) # since Poisson
        # Poisson uncertainty on observed rate
        #sigma_R = self.R_obs / np.sqrt(self.N_events)   # (N,)
        chi_R = -0.5 * np.sum((R_model_at_ztrue - self.R_obs)**2 / sigma_R**2)
        # --- Individual event log-likelihoods ---
        # L_i = integral R(z) * N(DL_obs_i - DL(z), sigma_DL_i) dz
        # Interpolate model DL onto catalogue z_true locations

        # Now both are at same locations — chi-squared is valid
        delta_DL = self.DL_obs - DL_model_at_ztrue 
        #delta_DL = self.DL_obs[:, None] - DL_grid[None, :]    # (N, Nz)
        sigma2   = sigma_DL[:, None]**2                        # (N, 1)
        chi_DL = -0.5 * np.sum(delta_DL**2/ sigma2)
        ''' 
        # --- DL vs z ---
        plt.scatter(self.z_true, np.log10(self.DL_obs), s=5, label='catalogue events')
        plt.plot(z_grid, np.log10(DL_grid), label='model DL(z)')
        plt.xlabel('z')
        plt.ylabel('DL [Mpc]')
        plt.show()
        
        p_z = np.maximum(R_grid, 0.0)
        p_z /= simpson(p_z, z_grid)
        plt.hist(self.z_true, bins=10, density=True, label='catalogue')
        plt.plot(z_grid, p_z, label='model p(z)')
        plt.xlabel('z')
        plt.ylabel('p(z)')
        plt.legend()
        plt.show()
        

        ln_gauss = -0.5 * np.sum(delta_DL**2/ sigma2)                # (N, Nz)
        #ln_R     = np.log(np.maximum(R_grid, 1e-300))         # (Nz,)
        
        ln_intgd = ln_gauss #+ ln_R[None, :]                # (N, Nz)

        dz       = z_grid[1] - z_grid[0]
        ln_max   = ln_intgd.max(axis=1, keepdims=True)
        ln_Li    = ln_max[:, 0] + np.log(
                       np.sum(np.exp(ln_intgd - ln_max) * dz, axis=1)
                   )   # (N,)
        '''
        if np.any(~np.isfinite(chi_DL + chi_R)):
            return -np.inf

        return  chi_DL + chi_R #-zeta + np.sum(ln_Li), np.sum(ln_R_tilde_ztrue)
