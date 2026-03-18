import numpy as np
import matplotlib.pyplot as plt


import reion_f as f
from reion_f import run_model

zstart, zend, dz = 30.0, 0.0, 0.2
Z_arraySize = int(round(abs((zend - zstart) / dz)))


(Z, QH_Q, dNLLdz, gamma_PI, 
         sfr_popII, sfr_popIII,
         dvc_dz, D_L, age_Gyr,
         tau_factor, ierr
            ) = run_model(
                h0     = 6.811000e+01,
                ombh2  = 2.260000e-02,
                omch2  = 1.179535e-01,
                ns     = 9.600000e-01,
                sigma_8= 0.8159,
                e_sf_ii= 0.01,
                alpha_z= 1e-05,
                esc_popii = 0.36,
                lambda0= 2.44,
                zstart_in=zstart,
                zend_in=zend,
                dz_in=dz,
                z_arraysize = Z_arraySize
            )

if ierr != 0:
    raise RuntimeError("filling() failed")


rho_sfrd_data = np.loadtxt('./ObsData/rho_SFRD_data.dat')

plt.plot(Z, QH_Q)
plt.xlim(2, 15)
plt.show()

plt.plot(Z, dNLLdz)
plt.xlim(2, 15)
plt.ylim(0, 15)
plt.show()

plt.plot(Z, gamma_PI*1e12)
plt.xlim(2, 15)
plt.ylim(0.01, 10.0)
plt.yscale('log')
plt.show()

plt.plot(Z, np.log10(sfr_popII), label='PopII')
#plt.plot(Z, np.log10(sfr_popIII), label='PopIII')
plt.scatter(rho_sfrd_data[:, 0], rho_sfrd_data[:, 1])
plt.xlim(5, 20)
plt.ylim(-4.5, -1.0)
plt.legend()
plt.show()

