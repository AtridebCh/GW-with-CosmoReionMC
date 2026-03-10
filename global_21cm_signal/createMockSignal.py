import sys

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

from brightnessTemp import Generate21cmSignal

YR_TO_SEC    = 3.1557e+07       # 

np.random.seed(42)

# ---------------------------------------------------------------------------
# Reionization model — Fortran f2py compiled module
# change the path
# ---------------------------------------------------------------------------
sys.path.append('../ReionizationWithF2Py/python_dir')
import reion_f as f
from reion_f import run_model

# ---------------------------------------------------------------------------
# Redshift grid — must match Fortran exactly
# ---------------------------------------------------------------------------
zstart, zend, dz = 30.0, 0.0, 0.2
Z_arraySize = int(round(abs((zend - zstart) / dz)))


# ---------------------------------------------------------------------------
# Run reionization model once at fiducial parameters
# ---------------------------------------------------------------------------

(Z_reion, QH_Q, dNLLdz, gamma_PI, 
         sfr_popII, sfr_popIII,
         dvc_dz, D_L, age_Gyr, 
         tau_factor, ierr) = run_model(
                h0=6.711000e+01,
                ombh2=2.260000e-02,
                omch2=1.179535e-01,
                ns=9.600000e-01,
                sigma_8=0.8022147706311814,
                esc_popii=0.0036,
                esc_popiii= 0.0,
                lambda0= 2.44,
                zstart_in=zstart,
                zend_in=zend,
                dz_in=dz,
                z_arraysize = Z_arraySize)

n_alpha_per_solar_mass = 9690*(1.12*1e57) #Furlanetto 06 provides N_alpha=9690 per unit baryon, multiplied by (M_Sun/M_p)
const_epsilon_x = 3.4*1e40 #chatterjee+ 19 in CGS

f_X     = 0.2
f_alpha = 1.0

Z_21cm = np.linspace(25.0, 6.0, 1000)


#sfr_popII in Msun/yr/Mpc^3

logEps_X     = np.log10(np.maximum(const_epsilon_x*sfr_popII, 1e-03)) #in erg/sec/mpc^3
lognDotAlpha = np.log10(np.maximum(n_alpha_per_solar_mass*sfr_popII/YR_TO_SEC, 1e-03)) #/sec/mpc^3

signal_setup = Generate21cmSignal(6.711000e+01, 2.260000e-02, 1.179535e-01, Z_reion,
                 logEps_X, lognDotAlpha, QH_Q, Z_21cm, f_X, f_alpha)
T_b          = signal_setup.signal_generator()

T_b += np.random.normal(0, 10.0, size=T_b.shape)   # adding 10 mK noise

plt.plot(Z_21cm, T_b)
plt.show()


np.savetxt('brightness_temp.txt', np.column_stack([Z_21cm, T_b]), 
           header='redshift T_b[mK]')
