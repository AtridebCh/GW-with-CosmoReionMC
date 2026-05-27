import numpy as np
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
from colossus.lss import mass_function

import reion_f as f
from reion_f import run_model, run_dndm, get_sfe

zstart, zend, dz = 30.0, 0.0, 0.2
Z_arraySize = int(round(abs((zend - zstart) / dz)))

cosmo = cosmology.setCosmology('planck18')

fzero    = 0.006, #0.17 
alpha_lo = 0.0, #0.65
alpha_hi = 0.0, #1.1

# Chatterjee. 2026; fzero    = 0.16, alpha_lo = 0.65, alpha_hi = 0.9, esc_popii = 0.2; works in new method
(Z, QH_Q, dNLLdz, gamma_PI, 
         sfr_popII, sfr_popIII,
         dvc_dz, D_L, age_Gyr,
         tau_factor, omega_dyn, omega_de, ierr
            ) = run_model(
                h0     = 6.811000e+01,
                ombh2  = 2.260000e-02,
                omch2  = 1.179535e-01,
                ns     = 9.600000e-01,
                sigma_8= 0.8159,
                omega_zero = -1.0,
                omega_a   = 0.0,
                fzero    = fzero, 
                alpha_lo = alpha_lo, 
                alpha_hi = alpha_hi,
                alpha_z  = 0.0,
                esc_popii = 0.25, #0.30
                lambda0= 2.64,
                zstart_in=zstart,
                zend_in=zend,
                dz_in=dz,
                z_arraysize = Z_arraySize
            )

if ierr != 0:
    raise RuntimeError("filling() failed")
    
'''
# Mass array and single redshift
m_arr = np.logspace(11, 15, 200)   # solar masses
z_val = np.array([7.0])            # single redshift as 1-element array

dndm_out, ierr = run_dndm(m_arr, z_val)



mfunc = mass_function.massFunction(m_arr*(cosmo.H0/100), 7.0, model = 'press74', q_out = 'dndlnM') #NOTE: colossus uses M-sun/h unit, so multiply by h for the conversion

plt.figure()
plt.loglog(m_arr, m_arr*np.abs(dndm_out[:, 0]), label ='cosmoreion')
plt.loglog(m_arr, mfunc*(cosmo.H0/100)**3, label = 'colossus')
plt.xlabel(r"$M\ [M_\odot]$")
plt.ylim(1e-07, 1e-01)
plt.ylabel(r"$dn/dlnM\ [\mathrm{Mpc}^{-3}]$")
plt.title(f"Halo Mass Function at z = {z_val[0]}")
plt.tight_layout()
plt.legend()
plt.savefig("dndm.pdf")
plt.show()

plt.plot(m_arr, m_arr*np.abs(dndm_out[:, 0])/(mfunc*(cosmo.H0/100)**3))
plt.show()


m_arr = np.logspace(11, 15.5, 200)   # solar masses
sfe_out, ierr = get_sfe(m_arr, fzero, alpha_lo, alpha_hi)

fig, ax = plt.subplots(figsize=(8, 6))


ax.plot(np.log10(m_arr), sfe_out)
ax.set_ylabel(r'$f_{*}(M_h)$')
ax.set_xlabel(r'$\log_{10} (M_h[M_{\odot}])$')
plt.tight_layout()
plt.savefig('f_star.pdf')
plt.show()


fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(Z, omega_de)
ax.set_xlim(2, 20)
ax.set_ylabel(r'$\Omega_{\rm DE}$')
ax.set_xlabel('redshift (z)')
plt.show()

'''

Lymanlimitdatafile = './ObsData/Lyman_limit.dat'
gammadatafile      = './ObsData/gamma_data_all_combined.dat'
# Ionizing background data
redshift_gamma, gamma, gamma_max, gamma_min = np.loadtxt(
            gammadatafile, usecols=(0, 1, 2, 3), unpack=True
        )
error_up       = gamma_max - gamma
error_down     = gamma- gamma_min

# Lyman-limit system data
lyman_redshift, lyman_limit_data, lyman_error = np.loadtxt(
            Lymanlimitdatafile, usecols=(0, 1, 2), unpack=True
        )


z_sfrd, sfrd, sfrd_plus, sfrd_minus = np.loadtxt(
    './ObsData/rho_SFRD_data.dat', usecols=(0, 1, 2, 3), unpack=True
)

# Asymmetric errors already in log space
error_up_sfrd   = sfrd_plus  - sfrd   # upper bar
error_down_sfrd = sfrd       - sfrd_minus  # lower bar


fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)

# --- Top left: QH_Q ---
ax = axes[0, 0]
ax.plot(Z, QH_Q, label='model', color='b')
#ax.set_xlim(2, 20)
ax.set_ylabel(r'$Q_{H}$')
ax.set_xlabel('redshift (z)')
#ax.text(12, 0.5, 'New Method')

# --- Top right: Lyman Limit ---
ax = axes[0, 1]
ax.errorbar(lyman_redshift, lyman_limit_data, yerr=lyman_error,
            fmt='o', capsize=3, label='data', color='k')

ax.plot(Z, dNLLdz, label='model', color='b')
#ax.set_xlim(2, 20)
ax.set_ylim(0.0, 12)
ax.set_ylabel(r'$dN_{LL}/dz$')
ax.set_xlabel('redshift (z)')
ax.legend()

# --- Bottom left: Gamma_PI ---
ax = axes[1, 0]
ax.errorbar(redshift_gamma, gamma, 
            yerr=[error_down, error_up ],   # asymmetric errors
            fmt='o', capsize=3, label='data', color='k')
ax.plot(Z, gamma_PI*1e12, label='model', color='b')
ax.set_yscale('log')
ax.set_ylim(1e-02, 1e01)
ax.set_ylabel(r'$\Gamma_{PI}$')
ax.set_xlabel('redshift (z)')
ax.legend()


# --- Bottom right: SFR ---
ax = axes[1, 1]
ax.errorbar(z_sfrd, sfrd,
            yerr=[error_down_sfrd, error_up_sfrd],
            fmt='o', capsize=3, label='data', color='k')
ax.plot(Z, np.log10(sfr_popII), label='PopII', color='b')
ax.set_ylim(-4.5, -1.0)
ax.set_ylabel(r'$\log_{10}(\rho_{SFR})$')
ax.set_xlabel('redshift (z)')
ax.legend()

plt.tight_layout()
plt.savefig('reionization_summary.pdf')
plt.show()
