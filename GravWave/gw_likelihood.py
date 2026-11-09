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
from scipy.stats import binned_statistic
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from project_paths import CATALOGUE_DIR

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
        catalogue_path = CATALOGUE_DIR / "CE_mock_catalog_13274_events_without_snr_noisy_sigmaD_x0.25_seed65.pkl", #"/Users/users/achatterjee/Data/GW_with_CosmoReionMC/GravWave/CE_mock_catalog_6636_events_with_snr_8_test.pkl", #  "/Users/users/achatterjee/Data/GW_with_CosmoReionMC/GravWave/CE_mock_catalog_6591_events_without_snr.pkl",
        #catalogue_path = "/Users/users/achatterjee/Data/GW_with_CosmoReionMC/GravWave/CE_mock_catalog_100_events_with_snr_8.pkl",
        pdet_path = CATALOGUE_DIR / "selection_function_p_det.pkl",
    ):
        self.catalogue_path = catalogue_path
        self.p_det_path = pdet_path
        self._is_setup      = False

    # ------------------------------------------------------------------
    # Setup — load catalogue once
    # ------------------------------------------------------------------

    def setup(self):
        """
        Load GW catalogue and precompute fiducial sigma_DL rescaling.

        Parameters
        ----------
        fiducial_H0, fiducial_Om0 : float
            Cosmology used when generating the catalogue.
        """
        with open(self.catalogue_path, "rb") as f:
            cat = pickle.load(f)
        '''
        self.N_events     = catalogue["observed_event"]
        print('N_events =', self.N_events)
        self.zeta_mock    = catalogue["expected_event"]
        #print('No. of events', self.N_events)
        self.DL_obs       = np.array([
            ev["injection"]["luminosity_distance"] for ev in catalogue
        ])                                                     # Mpc, shape (N,)
        self.sigma_DL_fid = np.array([ev["sigma_DL"] for ev in catalogue])
        self.z_true       = np.array([ev["z_true"] for ev in catalogue])
        self.R_obs        = np.array([ev["R_true"] for ev in catalogue])
        '''
        self.N_events     = cat["observed_event"]
        self.zeta_mock    = cat["expected_event"]
        self.T_obs        = cat["T_obs"]
        self.z_true      = np.array(cat["z_true"])
        self.R_obs      = np.array(cat["R_true"])


        a, b = cat["sigma_model"]["a"], cat["sigma_model"]["b"]
        self.sigma_of_D = lambda D: a * np.abs(D)**b
        #self.sigma_DL_fid    = np.array(cat["sigma_DL"])
        #injection   = cat["injection"]  # list of dicts, one per event
        #self.DL_obs = np.array([
        #                    inj["luminosity_distance"] for inj in cat["injection"]]) 
        self.DL_obs = np.array(cat["DL_obs"])
        
        #self.DL_obs_cat = np.array([inj["luminosity_distance"] for inj in cat["injection"]])
        #rng       = np.random.default_rng(7)
        #sigma_shf = rng.permutation(np.array(cat["sigma_DL"]))
        #self.sigma_DL_fid = sigma_shf
        #self.DL_obs = self.DL_obs_cat + rng.normal(0.0, sigma_shf)
        with open(self.p_det_path, "rb") as fpick:
             saved = pickle.load(fpick)

        self.p_det_interp = saved["p_det_interp"]

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

        #self.DL_obs = self.DL_true
        #np.savetxt('z_vs_R_from_mock.txt', np.c_[self.z_true, np.interp(self.z_true, z_grid, R_grid)])
        
        if R_grid is None or DL_grid is None:
            return -np.inf

        if np.any(~np.isfinite(R_grid)) or np.any(R_grid < 0):
            return -np.inf
        #R_model_at_z_true = np.interp(self.z_true, z_grid, R_grid)
        
        #norm = simpson(R_grid*self.p_det_interp(z_grid), z_grid)
        norm = simpson(R_grid, z_grid)
        zeta_model = self.T_obs * norm #should we multiply with 4*pi
        #sigma_zeta = 0.1*self.N_events
        #print('norm in zeta', norm)
        if zeta_model <= 0 or not np.isfinite(zeta_model):
            return -np.inf
         
        # --- Rescale sigma_DL for current cosmology ---
        DL_model_at_ztrue = np.interp(self.z_true, z_grid, DL_grid)
        #DL_model_atz6 = np.interp(6.0, z_grid, DL_grid)
        #sigma_DL          = self.sigma_DL_fid #* (
        #     DL_model_at_ztrue / self.DL_obs
        #) 
        # --- Individual event log-likelihoods ---
        # L_i = integral R(z) * N(DL_obs_i - DL(z), sigma_DL_i) dz
        # Interpolate model DL onto catalogue z_true locations

        # Now both are at same locations — chi-squared is valid
        #delta_DL = self.DL_obs[:, None] - DL_model_at_ztrue[:, None] #(N, 1)
        #sigma2   = sigma_DL**2
        # when we want to use ln_gauss
        
        #with noisy D_L
        #rng = np.random.default_rng(1) # seed = 1,2,3,4

        #self.DL_obs = self.DL_obs + rng.normal(0.0, self.sigma_DL_fid)
        #n_below = np.sum(self.DL_obs < DL_model_atz6)
        #print('n_below', n_below)
        

        sig_z    = self.sigma_of_D(DL_grid) 
        #delta_DL = self.DL_obs[:, None] - DL_grid[None, :]    # (N, Nz)
        
        #if you use sigma_DL = self.sigma_DL_fid *DL_model_at_ztrue / self.DL_obs) then
        # you have to add the ln_norm term in the ln_Li; else use sigma2   = self.sigma_DL_fid[:, None]**2
        #then you don't have to use the ln_norm in ln_Li
        #sigma2   = sigma_DL[:, None]**2 #self.sigma_DL_fid[:, None]**2                        # (N, 1)
        #ln_DL = -0.5*np.sum(delta_DL**2/ sigma2)
        ln_gauss = (-0.5 * ((self.DL_obs[:, None] - DL_grid[None, :]) / sig_z[None, :])**2
                    - np.log(sig_z)[None, :])
        #ln_gauss = -0.5*(delta_DL**2/ sigma2)                # (N, Nz)
        ln_R     = np.log(np.maximum(R_grid, 1e-300))         # (Nz,)
        #print(ln_R[0]+ ln_gauss[0,:])
        ln_intgd = ln_gauss + ln_R[None, :]  #(N, Nz)
        dz       = z_grid[1] - z_grid[0]
        ln_max   = ln_intgd.max(axis=1, keepdims=True)
        #ln_norm = -0.5*np.log(2*np.pi*sigma_DL**2) 
        integrand = np.exp(ln_intgd - ln_max)                              # (N, Nz)
        integral = simpson(integrand, z_grid, axis=1)
        ln_Li = ln_max[:, 0] + np.log(integral)
        #ln_Li = ln_max[:, 0] + np.log(integral) - np.log(norm) #for gaussian
         
        '''


        #check accuarcy
        difference = np.interp(self.z_true, z_grid, DL_grid) / self.DL_obs - 1
        print(np.min(difference), np.max(difference))
        '''
        '''
        # --- per-event integrals, each on its own grid starting at the edge ---
        sig   = self.sigma_DL_fid
        DL6   = DL_grid[0]                                         # D_L at z_grid[0] = 6
        u_e   = (DL6 - self.DL_obs) / sig                          # edge, in units of sigma
        u_lo  = np.maximum(u_e, -8.0)
        u_hi  = np.where(u_lo < 0, 8.0,
                         u_lo + np.minimum(16.0, 20.0 / np.maximum(u_lo, 1.0)))
        s     = np.linspace(0.0, 1.0, 201)
        u     = u_lo[:, None] + (u_hi - u_lo)[:, None] * s[None, :]    # (N, 201)

        DL_loc = self.DL_obs[:, None] + u * sig[:, None]
        valid  = DL_loc <= DL_grid[-1]
        z_loc  = np.interp(DL_loc, DL_grid, z_grid)
        J      = 1.0 / np.interp(z_loc, z_grid, np.gradient(DL_grid, z_grid))   # dz/dD_L
        R_loc  = np.interp(z_loc, z_grid, R_grid)

        ln_int = -0.5 * u**2 + np.log(np.maximum(R_loc, 1e-300)) + np.log(J)
        ln_int = np.where(valid, ln_int, -np.inf)
        m      = ln_int.max(axis=1, keepdims=True)
        ln_Li  = m[:, 0] + np.log(simpson(np.exp(ln_int - m), x=u, axis=1))
        '''
        '''
        '''

        '''
        #grid issue
        dDL_dz  = np.gradient(DL_grid, z_grid)
        #sigma_z = self.sigma_DL_fid / np.interp(self.z_true, z_grid, dDL_dz)
        sigma_z = self.sigma_of_D(self.DL_obs) / np.interp(self.z_true, z_grid, dDL_dz)
        #dz      = np.diff(z_grid)
        ratio   = sigma_z / (z_grid[1] - z_grid[0])
        print("min ratio:", ratio.min(), "  events below 1.3:", np.sum(ratio < 1.3), "of", len(ratio))
        #print("min sigma_z / dz =", (sigma_z / np.interp(self.z_true, z_grid[1:], dz)).min())
        '''
        
        '''
        #edge issue
        lo = (self.DL_obs - DL_grid[0])  / self.sigma_DL_fid    # sigmas above the lower edge
        hi = (DL_grid[-1] - self.DL_obs) / self.sigma_DL_fid    # sigmas below the upper edge
        print("closest to lower edge:", lo.min(), " events within 5 sigma:", np.sum(lo < 5))
        print("closest to upper edge:", hi.min(), " events within 5 sigma:", np.sum(hi < 5))
        print("mock z range :", self.z_true.min(), self.z_true.max())
        print("z_grid range :", z_grid[0], z_grid[-1])
        print("R at z_max / R peak:", R_grid[-1] / R_grid.max()) 
        '''
        #ln_Li = ln_max[:, 0] + np.log(integral) - np.log(norm) #for gaussian
        
        #count term
        #print(-zeta_model + self.N_events * np.log(zeta_model))
        
        #integral ratio
        #ratio   = np.interp(self.z_true, z_grid, DL_grid) / self.DL_true - 1
        #print("max |ratio|:", np.abs(ratio).max(), "  mean:", ratio.mean())
        
        #whether sigma_DL is added properly
        #r = (self.DL_obs - self.DL_true) / self.sigma_DL_fid
        #print("mean:", r.mean(), "  std:", r.std())  
        '''
        #Noise free D_L analysis
        z_i   = np.interp(self.DL_obs, DL_grid, z_grid)
        dDLdz = np.interp(z_i, z_grid, np.gradient(DL_grid, z_grid))
        R_i   = np.interp(z_i, z_grid, R_grid)
        count = -zeta_model + self.N_events * np.log(zeta_model)
        shape = np.sum(np.log(np.interp(z_i, z_grid, R_grid) / norm))
        jac   = -np.sum(np.log(np.interp(z_i, z_grid, np.gradient(DL_grid, z_grid))))
        print('count, shape, jac', count, shape, jac)
        '''
        if np.any(~np.isfinite(ln_Li)):
            return -np.inf
        #p = R_grid / simpson(R_grid, x=z_grid)
        #print(np.sum(np.log(np.interp(self.z_true, z_grid, p))))
        #return -zeta_model + np.sum(np.log(np.maximum(R_i, 1e-300)) - np.log(dDLdz))
        return -zeta_model + np.sum(ln_Li) #+np.log(self.T_obs) -self.DL_obs.shape[0]*np.log(zeta) #+ ln_DL
        #return ln_count + np.sum(ln_Li)
