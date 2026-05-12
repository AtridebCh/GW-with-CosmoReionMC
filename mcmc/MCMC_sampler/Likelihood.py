import sys
import os
import numpy as np
import logging
import matplotlib.pyplot as plt
from scipy import interpolate
from copy import deepcopy
import clik

from CoredataGenerator import CoreModule
from ChainContext import ChainContext
from utils import getLogger

# GW likelihood path, change this path
sys.path.append('../../GravWave')
from gw_likelihood import GWLikelihood


_CLIK_CACHE = {}

def _get_clik(path):
    """Load clik object once per process, cache by path."""
    if path not in _CLIK_CACHE:
        if path.endswith('.clik_lensing'):
            _CLIK_CACHE[path] = clik.clik_lensing(path)
        else:
            _CLIK_CACHE[path] = clik.clik(path)
    return _CLIK_CACHE[path]

class Likelihood:
    """
    Likelihood class for the reionization parameter estimation chain.

    Parameters
    ----------
    CoreModule  : callable
        Core module that computes model products.
    Planck_dict : dict
        Planck cosmological parameters.
    Astro_dict  : dict
        Astrophysical parameters.
    min_param   : array-like, optional
        Lower bounds for the free parameters.
    max_param   : array-like, optional
        Upper bounds for the free parameters.
    """

    def __init__(self, CoreModule,
                 min_param=None, max_param=None,
                 include_gw=True,
                 include_21cm = True):
        self.min        = min_param
        self.max        = max_param
        self.include_gw = include_gw
        self.include_21cm = include_21cm

    # ------------------------------------------------------------------
    # Parameter validation
    # ------------------------------------------------------------------

    def isValid(self, p):
        """Check whether parameter vector p lies within bounds."""
        if self.min is not None:
            for i in range(len(p)):
                if p[i] < self.min[i]:
                    getLogger().debug(
                        "Params out of bounds i=" + str(i) + " params " + str(p)
                    )
                    return False

        if self.max is not None:
            for i in range(len(p)):
                if p[i] > self.max[i]:
                    getLogger().debug(
                        "Params out of bounds i=" + str(i) + " params " + str(p)
                    )
                    return False

        return True

    # ------------------------------------------------------------------
    # Likelihood computation
    # ------------------------------------------------------------------

    def ComputeLikelihood(self, ctx, single_run=False):
        r"""
        Compute the total log-likelihood from model outputs stored in ctx.

        Parameters
        ----------
        ctx         : ChainContext
        single_run  : bool, optional

        Returns
        -------
        loglike_tot       : float
        tau               : float
        Q_HII_at_z5p8     : float
        """
        redshift         = ctx.get('Redshift')
        Q_HII            = ctx.get('QHII_array')
        tau              = ctx.get('tau')
        Q_HII_at_z5p8    = ctx.get('Q_HII_at_Z5.8')
        lyman_limit      = ctx.get('LymanLimit')
        gamma_PI         = ctx.get('Gamma_PI')
        cl_TT            = ctx.get('cl_TT')
        cl_EE            = ctx.get('cl_EE')
        cl_BB            = ctx.get('cl_BB')
        cl_TE            = ctx.get('cl_TE')
        

        ells = np.arange(len(cl_TT)) + 2
        cls = np.array([2.*np.pi*cl_TT/(ells*(ells+1.)), 2.*np.pi*cl_EE/(ells*(ells+1.)),
                     2.*np.pi*cl_BB/(ells*(ells+1.)), 2.*np.pi*cl_TE/(ells*(ells+1.)),
                     None, None], dtype='object')
        p = deepcopy(ctx.getParams())
        cl01 = np.array([0.0, 0.0])

        if 0.135>tau>0.02:# and Q_HII_at_z5p8>0.94:  #dark pixel fraction gives an upper limit of x_HI=0.06 hence QHII>(1-0.06)
            loglike = 0.

            for i, clik_file in enumerate(self.clik_files):
                l = _get_clik(clik_file)    # loaded once per worker, cached after that
                inds = np.array([bool(int(f)) for f in l.has_cl])
                lmax = np.array(l.lmax)
                input = np.array([])
                
                for j, flag in enumerate(inds):
                    if flag:
                        tcl = np.append(cl01, cls[j][:lmax[j]-1])
                        input = np.append(input, tcl)
                if self.clik_files[i]=='../../COM_Likelihood_Data-baseline_R3.00/baseline/plc_3.0/low_l/commander/commander_dx12_v3_2_29.clik':
                    nuisense=[1.0]
                if self.clik_files[i]=='../../COM_Likelihood_Data-baseline_R3.00/baseline/plc_3.0/hi_l/plik_lite/plik_lite_v22_TTTEEE.clik':
                    nuisense=[1.0]
                if self.clik_files[i]=='../../COM_Likelihood_Data-baseline_R3.00/baseline/plc_3.0/low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik':                   
                    nuisense=[1.0] 
                input = np.hstack([input, nuisense])
                loglike += l(input)[0]
            
            # redshift is descending (30→0) — must reverse for np.interp
            lyman_limit_interpolated = np.interp(
                                       self.lyman_redshift,
                                       redshift[::-1],      
                                       lyman_limit[::-1]   
                                       )

            gamma_interpolated = np.interp(
                                       self.redshift_gamma,
                                       redshift[::-1],          
                                       np.log10(gamma_PI)[::-1]
                                       )
            
            loglike_lyman  = -0.5 * np.sum(
                             (self.lyman_limit_data - lyman_limit_interpolated) ** 2 / self.lyman_error ** 2
                           )
            loglike_gamma  = -0.5 * np.sum(
                              (self.gamma_log - gamma_interpolated) ** 2 / self.error_log ** 2
                           )
            loglike_tau    = -0.5 * (tau - 0.054) ** 2 / (0.007 ** 2)
            loglike_tot    = loglike + loglike_lyman + loglike_gamma + loglike_tau
            
            # ---- GW hyperlikelihood ----
            if self.include_gw:
                loglike_gw   = self.gw_like.ComputeLikelihood(ctx)
                loglike_tot += loglike_gw
                getLogger().debug(f"ln_L_GW = {loglike_gw:.3f}")
            #----- 21cm liklihood
            if self.include_21cm:
                T_b_model = ctx.get('T_b')
                loglike_21cm = -0.5 * np.sum(
                               (self.T_b_Obs - T_b_model) ** 2 / self.sigma_T_b ** 2
                               )
                loglike_tot += loglike_21cm
            

            if single_run:
                print('loglike cmb, loglike_lyman, loglike_gamma, \
                      loglike_tau, loglike_gw, loglike_21cm, loglike_tot', loglike, loglike_lyman, 
                                        loglike_gamma, loglike_tau, loglike_gw, loglike_21cm, loglike_tot)
                print('tau=', tau)
                
                plt.plot(redshift[::-1], Q_HII[::-1])
                #plt.plot(self.redshift_gamma, gamma_interpolated)
                plt.xlabel('redshift (z)')
                plt.ylabel(r'$Q_{mathrm{HII}}$')
                plt.savefig('Q_HII.pdf')
                plt.show()

                 
                plt.scatter(self.redshift_gamma, self.gamma_log)
                plt.plot(self.redshift_gamma, gamma_interpolated)
                plt.xlabel('redshift (z)')
                plt.ylabel(r'$\Gamma_{PI}$')
                plt.savefig('Gamma_PI.pdf')
                plt.show()
                
                plt.scatter(self.lyman_redshift, self.lyman_limit_data)
                plt.plot(self.lyman_redshift, lyman_limit_interpolated)
                plt.xlabel('redshift (z)')
                plt.ylabel(r'$dN_{LL}/dz$')
                plt.savefig('Lyman_Limit.pdf')
                plt.show()
                
                plt.plot(self.T_b_redshift, self.T_b_Obs)
                plt.plot(self.T_b_redshift, T_b_model)
                plt.xlabel('redshift (z)')
                plt.ylabel(r'$T_{b}(mk)$')
                plt.savefig('T_b.pdf')
                plt.show()

            return loglike_tot, tau, Q_HII_at_z5p8
        else:
            return -np.inf, tau, Q_HII_at_z5p8

    # ------------------------------------------------------------------
    # Call interface for sampler
    # ------------------------------------------------------------------

    def __call__(self, p, single_run=False):
        getLogger().debug("pid: %s, processing: %s" % (os.getpid(), p))
        ctx = self.createChainContext(p)

        if np.all(p > 0.0):
            #calling Core to compute all observables for current mcmc step
            model = self.Core(ctx, single_run=single_run)
            if model:
                loglike, tau, Q_HII_at_z5p8 = self.ComputeLikelihood(
                    ctx, single_run=single_run
                )
                blobs = np.array([tau, Q_HII_at_z5p8])
                return loglike, blobs
            else:
                return -np.inf, [1.0, np.nan]
        else:
            return -np.inf, [1.0, np.nan]

    # ------------------------------------------------------------------
    # Chain context factory
    # ------------------------------------------------------------------

    def createChainContext(self, p):
        """Return a new ChainContext instance for parameter vector p."""
        try:
            p = Params(*zip(self.params.keys, p))
        except Exception:
            # no params object or params has no keys
            pass
        return ChainContext(self, p)

    # ------------------------------------------------------------------
    # Setup — load observational data and initialise core module
    # ------------------------------------------------------------------

    def setup(self, gamma_datafile, lyman_limit_datafile, T_b_Datafile, single_run=False):
        """
        Load observational data files and initialise the core module.

        Parameters
        ----------
        gamma_datafile       : str
            Path to ionizing background (Gamma) data file.
        lyman_limit_datafile : str
            Path to Lyman-limit system data file.
        """
        # Ionizing background data
        self.redshift_gamma, gamma, gamma_max, gamma_min = np.loadtxt(
            gamma_datafile, usecols=(0, 1, 2, 3), unpack=True
        )
        self.gamma_log = np.log10(gamma)
        error_up       = gamma_max - gamma
        self.error_log = error_up / gamma

        # Lyman-limit system data
        self.lyman_redshift, self.lyman_limit_data, self.lyman_error = np.loadtxt(
            lyman_limit_datafile, usecols=(0, 1, 2), unpack=True
        )
        
        #change these paths
        self.clik_files = [
            "../../COM_Likelihood_Data-baseline_R3.00/baseline/plc_3.0/low_l/commander/commander_dx12_v3_2_29.clik",
            "../../COM_Likelihood_Data-baseline_R3.00/baseline/plc_3.0/low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik",
            "../../COM_Likelihood_Data-baseline_R3.00/baseline/plc_3.0/hi_l/plik_lite/plik_lite_v22_TTTEEE.clik",
        ]

        self.Core = CoreModule()
        
        # GW likelihood — initialised once, catalogue loaded once
        if self.include_gw:
            self.gw_like = GWLikelihood()
            self.gw_like.setup()
        if self.include_21cm:
           self.T_b_redshift, self.T_b_Obs = np.loadtxt(T_b_Datafile, unpack =True)
           self.sigma_T_b = 10 #mk
