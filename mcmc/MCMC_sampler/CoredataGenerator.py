import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import simpson
import sys

from warnings import filterwarnings
filterwarnings('ignore')

import camb

# we have checked that importing outside is not creating any problem
sys.path.append('../../ReionizationWithF2Py/python_dir') #/home/atri/GW_with_CosmoReionMC/ReionizationWithF2Py/python_dir
import reion_f as f
from reion_f import run_model, run_dndm, get_sfe

sys.path.append('../../GravWave')
from gw_model import GWModel


sys.path.append('../../global_21cm_signal')
from brightnessTemp import Generate21cmSignal

YR_TO_SEC    = 3.1557e+07       # yr/s

# Free parameter name → index mapping
FREE_PARAM_MAPPING_DEFAULT = {
    'H0':      0,
    'ombh2':   1,
    'omch2':   2,
    'As':      3,
    'ns':      4,
    'omega_zero': 5,
    'omega_a': 6,
    'fzero':   7,
    'alpha_lo': 8,
    'alpha_hi': 9,
    'alpha_z': 10, 
    'esc_II' : 11,
    'lambda0': 12,
    'f_X'    : 13,
    'f_alpha': 14
}

# Redshift grid for reionization
zstart, zend, dz = 30.0, 0.0, 0.2
Z_arraySize = int(round(abs((zend - zstart) / dz)))

# exact match to Fortran as reionization code uses f2py
#Z = np.array([ik * (-abs(dz)) + zstart for ik in range(Z_arraySize + 1)])

#if you change this and t_delay, you have to change in the GWcatalogue and 21cm signal too 
#and create the catalogues agian
#so best left unchanged
Z_21cm = np.linspace(25.0, 6.0, 1000)
Z_GW   = np.linspace(6.0, 25.0, 500)
t_delay = 0.1

CAMB_KMIN = -4    # log10(k_min) in h/Mpc; Works fine wiht this value
CAMB_KMAX =  1    # log10(k_max) in h/Mpc; works fine with this value

n_alpha_per_solar_mass = 9690*(1.12*1e57) #Furlanetto 06 provides N_alpha=9690 per unit baryon, multiplied by (M_Sun/M_p)
const_epsilon_x = 3.4*1e40 #chatterjee+ 19 in CGS

# ---------------------------------------------------------------------------
# GWModel instantiated once at module level — not per likelihood call
# ---------------------------------------------------------------------------
_gw_model = GWModel(
    Z_GW,
    t_delay = t_delay,
) 

class CoreModule:
    """
    Core module for delegating computation of model products.

    Parameters
    ----------
    params_mapping : dict
        Maps free parameter names to their index in the parameter vector.
        Defaults to FREE_PARAM_MAPPING_DEFAULT.
    """

    def __init__(self, params_mapping=None):
        if params_mapping is None:
            params_mapping = FREE_PARAM_MAPPING_DEFAULT
        self.params_mapping = params_mapping
        self.lmax = 2508
        #self.reion_path     = '../../ReionizationWithF2Py/python_dir'

    def __call__(self, ctx, single_run=False):
        p = ctx.getParams()

        try:
            # ----------------------------------------------------------
            # Unpack free parameters
            # ----------------------------------------------------------
            free_params = {}
            for name, idx in self.params_mapping.items():
                free_params[name] = p[idx]

            # ----------------------------------------------------------
            # CAMB k-grid (h/Mpc)
            # ----------------------------------------------------------
            self._k_camb = np.logspace(CAMB_KMIN, CAMB_KMAX, 200)  # do not change

            # ----------------------------------------------------------
            # First CAMB call — base cosmology, no reionization
            # Purpose: compute sigma8 to pass into reionization model
            # ----------------------------------------------------------
            camb_params_base = camb.CAMBparams(
                Reion=camb.reionization.TanhReionization(Reionization=False)
            )
            camb_params_base.set_cosmology(
                H0=free_params['H0'],
                ombh2=free_params['ombh2'],
                omch2=free_params['omch2'],
            )
            camb_params_base.InitPower.set_params(
                ns=free_params['ns'],
                As=free_params['As'],
            )
            camb_params_base.set_matter_power(
                redshifts=[0.0],
                kmax=CAMB_KMAX,
            )
            camb_params_base.set_dark_energy(
                w=free_params['omega_zero'], 
                wa=free_params['omega_a'], 
                dark_energy_model='fluid',
            )

            camb_results_base = camb.get_results(camb_params_base)
            sigma8 = camb_results_base.get_sigma8()[-1]  # sigma8 at z=0
            

            # ----------------------------------------------------------
            # Run reionization model using sigma8 from first CAMB call
            
            (
                Z_reion,
                QH_Q, 
                dNLLdz, 
                gamma_PI, 
                sfr_popII, 
                sfr_popIII,
                dvc_dz, 
                D_L, 
                age_Gyr, 
                tau_factor,
                omega_dyn, 
                omega_de,
                ierr
            ) = run_model(
                h0=free_params['H0'],
                ombh2=free_params['ombh2'],
                omch2=free_params['omch2'],
                ns=free_params['ns'],
                sigma_8=sigma8,
                omega_zero = free_params['omega_zero'],
                omega_a = free_params['omega_a'],
                fzero =free_params['fzero'],
                alpha_lo = free_params['alpha_lo'], 
                alpha_hi = free_params['alpha_hi'],
                alpha_z= 0.0,
                esc_popii = free_params['esc_II'],
                lambda0=free_params['lambda0'],
                zstart_in=zstart,
                zend_in=zend,
                dz_in=dz,
                z_arraysize = Z_arraySize
            )
            
            idx_z5p8      = np.argmin(np.abs(Z_reion - 5.8))
            Q_HII_at_z5p8 = QH_Q[idx_z5p8]
            
            plt.plot(Z_reion, QH_Q)
            plt.show()


            # ----------------------------------------------------------
            # Second CAMB call — inject external reionization history
            # ----------------------------------------------------------
            camb_params_reion = camb.CAMBparams()
            camb_params_reion.set_cosmology(
                H0=free_params['H0'],
                ombh2=free_params['ombh2'],
                omch2=free_params['omch2'],
            )
            camb_params_reion.InitPower.set_params(
                ns=free_params['ns'],
                As=free_params['As'],
            )
            camb_params_reion.set_matter_power(
                redshifts=[0.0],
                kmax=self._k_camb[-1],
            )
            camb_params_reion.set_dark_energy(
                 w=free_params['omega_zero'], 
                 wa=free_params['omega_a'], 
                 dark_energy_model='fluid',
            )
            camb_params_reion.NonLinear          = camb.model.NonLinear_none
            camb_params_reion.ReionExternal      = True
            camb_params_reion.Reion.UsePCReion   = True
            camb_params_reion.length_array       = 151
            camb_params_reion.redshift_array_external = Z_reion[:-1]
            camb_params_reion.reionization_history    = QH_Q[:-1]
            camb_params_reion.set_for_lmax(self.lmax, lens_potential_accuracy=2)


            # ----------------------------------------------------------
            # Compute CMB power spectra
            # ----------------------------------------------------------
            results = camb.get_results(camb_params_reion)
            powers  = results.get_cmb_power_spectra(
                         lmax    = self.lmax,
                         CMB_unit = 'muK',
                         raw_cl  = False      # returns Dℓ = ℓ(ℓ+1)Cℓ/2π in μK²
                         )['total']

            data_bg= camb.get_background(camb_params_reion)
            #redshift_tau= np.linspace(0.0, 10.0, 101) #redshift for calculating tau, no need to go beyond 10
            
            back_ev = data_bg.get_background_redshift_evolution(Z_reion, ['x_e'])
            tau     = simpson(back_ev['x_e'] * tau_factor, Z_reion) 
            
            
                  
            R_merger, DL_grid, z_grid_gw = _gw_model.compute(free_params['H0'], free_params['ombh2'],
                                              free_params['omch2'], Z_reion, sfr_popII, age_Gyr, D_L, dvc_dz)
            
            logEps_X     = np.log10(np.maximum(const_epsilon_x*sfr_popII, 1e-03)) #in erg/sec/mpc^3
            lognDotAlpha = np.log10(np.maximum(n_alpha_per_solar_mass*sfr_popII/YR_TO_SEC, 1e-03)) #/sec/mpc^3                              
            signal_setup = Generate21cmSignal(free_params['H0'], free_params['ombh2'],
                                              free_params['omch2'], Z_reion, logEps_X, 
                                              lognDotAlpha, QH_Q, Z_21cm, free_params['f_X'], free_params['f_alpha'])
            T_b          = signal_setup.signal_generator()

            # ----------------------------------------------------------
            # camb and reionization outputs in ctx
            # ----------------------------------------------------------
            ctx.add('cl_TT',    powers[2:, 0])
            ctx.add('cl_EE',    powers[2:, 1])
            ctx.add('cl_BB',    powers[2:, 2])
            ctx.add('cl_TE',    powers[2:, 3])
            ctx.add('Redshift', Z_reion)
            ctx.add('QHII',     QH_Q)
            ctx.add('tau',      tau) 
            ctx.add('Q_HII_at_Z5.8',  Q_HII_at_z5p8)
            ctx.add('LymanLimit',      dNLLdz)
            ctx.add('Gamma_PI',        gamma_PI*1e12) #converting to gamma_PI_{-12} unit
            
            # --- GW in ctx ---
            ctx.add('gw_R_merger', R_merger)
            ctx.add('gw_DL_grid',  DL_grid)
            ctx.add('gw_z_grid',   z_grid_gw)
            
            #--- 21cm in ctx ---
            ctx.add('T_b', T_b)

            return 1.0

        except Exception:
            raise Exception(
                "Error either in assigning values in ctx or while "
                "computing some quantity from the model."
            )
